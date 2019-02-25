/*
(c) 2019 Christian Henke
Licensed under Apache License 2.0
https://www.apache.org/licenses/LICENSE-2.0


Convert fasta to Elasticsearch JSON while including contig-based annotation information and coverage.
*/

package main

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"encoding/csv"
	"encoding/json"
	"flag"
	"fmt"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/sam"
	"github.com/shenwei356/bio/seqio/fastx"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
	"time"
)

const CovAvgSpan int = 10

type Contig struct {
	ID         string     `json:"id"`
	Nucleotide string     `json:"nucleotide"`
	Length     int        `json:"length"`
	Coverage   []Coverage `json:"coverage,omitempty"`
	Virsorter  *Virsorter `json:"virsorter,omitempty"`
	Rpkm       []Rpkm     `json:"rpkm,omitempty"`
}

type Coverage struct {
	Name   string `json:"name"`
	Values []int  `json:"values"`
}

type Virsorter struct {
	Phage    string `json:"phage,omitempty"`
	Prophage string `json:"prophage,omitempty"`
}

type Rpkm struct {
	Name  string  `json:"name"`
	Value float64 `json:"value"`
}

type Contigs map[string]*Contig

var (
	fastaFile       = flag.String("fasta", "", "FASTA input file")
	virSorterCSV    = flag.String("virsorter-csv", "", "VirSorter input file")
	sampleNames     = flag.String("sample-names", "", "Comma-separated list of sample names")
	sampleBAMFiles  = flag.String("sample-bam-files", "", "Comma-separated list of sorted BAM files")
	sampleRPKMFiles = flag.String("sample-rpkm-files", "", "Comma-separated list of BBMap .rpkm files (gzipped)")
	jsonOutputFile  = flag.String("json", "", "JSON output file")
	skipCoverage    = flag.Bool("skip-coverage", false, "Skip coverage calculation")
	help            = flag.Bool("help", false, "Display help")
)

func main() {
	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}

	// check for required inputs

	if *fastaFile == "" {
		log.Fatal("Error: FASTA input file parameter is required.")
	}
	if *jsonOutputFile == "" {
		log.Fatal("Error: JSON output file parameter is required.")
	}
	var sampleNameList, bamFileNames, rpkmFileNames []string
	if *sampleBAMFiles != "" || *sampleRPKMFiles != "" {
		sampleNameList = strings.Split(*sampleNames, ",")
	}
	if *sampleBAMFiles != "" {
		if *sampleNames == "" {
		}
		bamFileNames = strings.Split(*sampleBAMFiles, ",")
		if len(sampleNameList) != len(bamFileNames) {
			log.Fatal("Error: Sample names / BAM files arguments differ in length.")
		}
	}
	if *sampleRPKMFiles != "" {
		if *sampleNames == "" {
			log.Fatal("Error: Sample names required for RPKM file parsing.")
		}
		rpkmFileNames = strings.Split(*sampleRPKMFiles, ",")
		if len(sampleNameList) != len(rpkmFileNames) {
			log.Fatal("Error: Sample names / RPKM files arguments differ in length.")
		}
	}

	start := time.Now()

	contigs := make(Contigs)

	// FASTA
	log.Print("Parsing contigs FASTA file ...")
	parseContigsFastaFile(*fastaFile, &contigs)

	// BAM
	if !*skipCoverage && *sampleBAMFiles != "" {
		log.Print("Parsing BAM files ...")
		for sampleIndex, sampleName := range sampleNameList {
			parseBAMFile(sampleName, bamFileNames[sampleIndex], &contigs)
		}
	}

	// RPKM
	if *sampleRPKMFiles != "" {
		log.Print("Parsing RPKM files ...")
		for sampleIndex, sampleName := range sampleNameList {
			parseRPKMFile(sampleName, rpkmFileNames[sampleIndex], &contigs)
		}
	}

	// VirSorter
	if *virSorterCSV != "" {
		log.Print("Parsing VirSorter file ...")
		parseVirSorterCSV(*virSorterCSV, &contigs)
	}

	// JSON
	log.Print("Writing JSON output file ...")
	objCount := writeJson(&contigs, *jsonOutputFile)
	log.Printf("%d JSON objects written to %q.", objCount, *jsonOutputFile)

	elapsed := time.Since(start)
	log.Printf("Total time: %s", elapsed)
}

func headerToID(header string) string {
	headerFields := strings.SplitN(header, " ", 2)
	if len(headerFields) < 2 || len(headerFields[0]) == 0 {
		log.Fatalf("Error while parsing header line: %q", header)
	}
	return headerFields[0]
}

func parseContigsFastaFile(path string, contigs *Contigs) {
	reader, err := fastx.NewDefaultReader(path)
	if err != nil {
		log.Fatalf("error reading FASTA: %v\n", err)
	}
	defer reader.Close()
	for {
		fastaRec, err := reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			if err != nil {
				log.Fatalf("error reading FASTA: %v\n", err)
			}
			break
		}
		(*contigs)[string(fastaRec.ID)] = &Contig{
			ID:         string(fastaRec.ID),
			Nucleotide: string(fastaRec.Seq.Seq),
			Length:     fastaRec.Seq.Length(),
		}
	}
}

func openBAMFile(path string) (*os.File, *bam.Reader) {
	var r io.Reader
	f, err := os.Open(path)
	if err != nil {
		log.Fatalf("Could not open file %q\n", err)
	}
	ok, err := bgzf.HasEOF(f)
	if err != nil {
		log.Fatalf("Could not open file %q\n", err)
	}
	if !ok {
		log.Printf("File %q has no BGZF magic block: may be truncated\n", path)
	}
	r = f

	b, err := bam.NewReader(r, 0)
	if err != nil {
		log.Fatalf("Could not read BAM: %q\n", err)
	}
	return f, b
}

func parseBAMFile(sampleName string, path string, contigs *Contigs) {
	f, b := openBAMFile(path)
	defer f.Close()
	defer b.Close()
	b.Omit(bam.AllVariableLengthData)

	var lastSeenContig *sam.Reference
	var cov []int
	for {
		bamRec, err := b.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("error reading BAM: %v\n", err)
		}
		if lastSeenContig == nil {
			lastSeenContig = bamRec.Ref
			cov = make([]int, bamRec.Ref.Len())
		}
		if bamRec.Ref != lastSeenContig {
			// wrap up current contig and save results
			contigID := headerToID(lastSeenContig.Name())
			c := Coverage{Name: sampleName, Values: calcCovAvg(&cov)}
			(*contigs)[contigID].Coverage = append((*contigs)[contigID].Coverage, c)

			// stop parser loop at the first unmapped sequence
			if bamRec.Ref.Name() == "*" {
				break
			}

			// prepare for mappings of next contig
			lastSeenContig = bamRec.Ref
			cov = make([]int, bamRec.Ref.Len())
		}
		if bamRec.Ref == lastSeenContig {
			pileupRec(bamRec, &cov)
		}
	}
}

func pileupRec(rec *sam.Record, cov *[]int) {
	recStart := rec.Start()
	recEnd := rec.End()
	for i := 0; i < rec.Ref.Len(); i++ {
		if recStart <= i && recEnd >= i {
			(*cov)[i]++
		}
	}
}

func calcCovAvg(cov *[]int) []int {
	covAvg := make([]int, 0, len(*cov)/CovAvgSpan+1)
	for i := 0; i < len(*cov)-CovAvgSpan; i += CovAvgSpan {
		tmp := 0
		for j := i; j < i+CovAvgSpan; j++ {
			tmp += (*cov)[j]
		}
		covAvg = append(covAvg, tmp/CovAvgSpan)
	}
	return covAvg
}

func parseRPKMFile(sampleName string, path string, contigs *Contigs) {
	tsvGz, err := os.Open(path)
	if err != nil {
		log.Fatalf("Could not open file: %q\n", err)
	}
	defer tsvGz.Close()
	tsv, err := gzip.NewReader(tsvGz)
	if err != nil {
		log.Fatalf("Could not open gzipped file: %q\n", err)
	}
	defer tsv.Close()

	// skip header rows
	scanner := bufio.NewScanner(tsv)
	scanner.Scan()
	scanner.Scan()
	scanner.Scan()
	scanner.Scan()
	scanner.Scan()
	tsvContents := new(bytes.Buffer)
	for scanner.Scan() {
		tsvContents.WriteString(scanner.Text())
		tsvContents.WriteString("\n")
	}

	tsvReader := csv.NewReader(tsvContents)
	tsvReader.Comma = '\t'
	for {
		record, err := tsvReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("Error parsing RPKM file: %q\n", err)
		}
		contigID := headerToID(record[0])
		val, err := strconv.ParseFloat(record[5], 64)
		if err != nil {
			log.Fatalf("Error parsing RPKM value: %q\n", err)
		}
		r := Rpkm{Name: sampleName, Value: val}
		(*contigs)[contigID].Rpkm = append((*contigs)[contigID].Rpkm, r)
	}
}

func parseVirSorterCSV(path string, contigs *Contigs) {
	vsFile, err := os.Open(path)
	if err != nil {
		log.Fatalf("Could not open file %q: %q\n", path, err)
	}
	defer vsFile.Close()
	var phageConfidence, prophageConfidence, state string
	scanner := bufio.NewScanner(vsFile)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		switch line {
		case "## 1 - Complete phage contigs - category 1 (sure)":
			phageConfidence = "sure"
			state = "phage"
		case "## 2 - Complete phage contigs - category 2 (somewhat sure)":
			phageConfidence = "somewhat sure"
			state = "phage"
		case "## 3 - Complete phage contigs - category 3 (not so sure)":
			phageConfidence = "not so sure"
			state = "phage"
		case "## 4 - Prophages - category 1 (sure)":
			prophageConfidence = "sure"
			state = "prophage"
		case "## 5 - Prophages - category 2 (somewhat sure)":
			prophageConfidence = "somewhat sure"
			state = "prophage"
		case "## 6 - Prophages - category 3 (not so sure)":
			prophageConfidence = "not so sure"
			state = "prophage"
		default:
			if !strings.HasPrefix(line, "#") {
				virsorter := &Virsorter{}
				contigID := strings.SplitN(strings.TrimPrefix(strings.SplitN(line, ",", 2)[0], "VIRSorter_"), "_flag=", 2)[0]
				if len(contigID) == 0 || strings.ContainsAny(contigID, "=") {
					log.Fatalf("VirSorter contig identifiers are being parsed incorrectly. "+
						"Original ID is %q, parsed contigID is %q", strings.SplitN(line, ",", 2)[0], contigID)
				}
				switch state {
				case "phage":
					virsorter.Phage = phageConfidence
				case "prophage":
					virsorter.Prophage = prophageConfidence
				}
				(*contigs)[contigID].Virsorter = virsorter
			}
		}

	}
}

func writeJson(contigs *Contigs, path string) int {
	out, err := os.Create(*jsonOutputFile)
	if err != nil {
		log.Fatalf("Could not open file for writing: %q\n", err)
	}
	defer out.Close()
	var i int
	for k, v := range *contigs {
		_, err := fmt.Fprintf(out, "{\"index\": {\"_id\": %s}}\n", k)
		if err != nil {
			log.Fatalf("Could not write to file: %q\n", err)
		}
		obj, err := json.Marshal(v)
		if err != nil {
			log.Fatalf("Could not generate JSON object for key %s:", k)
		}
		_, err = fmt.Fprintf(out, "%s\n", string(obj))
		if err != nil {
			log.Fatalf("Could not write to file: %q\n", err)
		}
		i++
	}
	return i
}

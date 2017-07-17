#!/usr/bin/env perl
#
# v1.0
# (c) Alexander Sczyrba
# Licensed under Apache License 2.0
# https://www.apache.org/licenses/LICENSE-2.0
#
# Modified by Christian Henke

use Getopt::Long;
use strict;
use Bio::SearchIO;
use Bio::LITE::Taxonomy::NCBI;

my $taxfile;
my $dbpath;

GetOptions("f=s"     => \$taxfile,
           "db=s"     => \$dbpath,
           );

unless ($taxfile) {
    die "\nusage: $0 -f TAXFILE -db DBPATH

reads a TAB-delimited file with GENEID TAXID and
translates the TAXID to complete lineage string
";
}
unless ($dbpath) {
    die "\nusage: $0 -f TAXFILE -db DBPATH

reads a TAB-delimited file with GENEID TAXID and
translates the TAXID to complete lineage string
";
}

print STDERR "reading Taxonomy DB...";
my $taxDB = Bio::LITE::Taxonomy::NCBI->new (
    db=>"NCBI",
    names=> "$dbpath/names.dmp",
    nodes=>"$dbpath/nodes.dmp",
#    dict=>"$dbpath/gi_taxid_prot.bin"
    );
print STDERR "done.\n";

open(IN,$taxfile) || die "Cannot open TAXFILE: $!";
while(<IN>) {
    chomp;
    my $geneid;
    my $taxid;
    ($geneid,$taxid) = split /\t/;
#    next if ($taxid == 1);
#    next if ($taxid == 2);
#    next if ($taxid == 131567);

    print "$geneid\t$taxid";

    my $bestlevel = "unknown";
    my $besttaxon = "unknown";
    my $lineage = "";

    foreach my $level ("superkingdom","phylum","class","order","family","genus","species","subspecies") {
        my $taxon = $taxDB->get_term_at_level($taxid,$level);
        if (($taxon ne "undef") && ($taxon ne "")) {
            $lineage .= "\t".$taxon;
            $bestlevel = $level;
            $besttaxon = $taxon;
        }
        else {
            $lineage .= "\tuc_".$besttaxon;
        }
    }
    print "\t$bestlevel\t$besttaxon\t$lineage\n";
}

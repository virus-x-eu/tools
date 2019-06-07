# Enter this directory
cd sra-metadata-search/

# Download metadata from https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/
wget -P /tmp/ "https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/NCBI_SRA_Metadata_Full_20190504.tar.gz"

# Convert metadata XML to JSON
./sra-metadata-to-json.py /tmp/NCBI_SRA_Metadata_Full_20190504.tar.gz data/submissions.gz

# Start Elasticsearch
docker-compose up

# Import submissions
docker-compose exec import /import/import.sh

# Open http://localhost
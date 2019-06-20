# Enter this directory
cd sra-metadata-search/

# Download metadata from https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/
wget -P /tmp/ "https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/NCBI_SRA_Metadata_Full_20190504.tar.gz"

# Create data directories
mkdir data
mkdir esdata
chmod 777 esdata

echo "vm.max_map_count=262144" | sudo tee /etc/sysctl.conf
sudo sysctl -f

# Convert metadata XML to JSON
./sra-metadata-to-json.py /tmp/NCBI_SRA_Metadata_Full_20190504.tar.gz data/submissions.gz

# Start Elasticsearch
docker-compose up

# Import submissions
docker-compose exec import /import/import.sh

# Open http://localhost

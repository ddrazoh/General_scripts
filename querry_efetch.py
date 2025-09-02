from Bio import Entrez
from tqdm import tqdm
import csv
import time
import xml.etree.ElementTree as ET

def fetch_shigella_genome_accessions_and_publications(email, db="nucleotide", term="Shigella[Organism] AND complete genome", retmax=2):
    try:
        Entrez.email = email

        # Search the NCBI nucleotide database for Shigella complete genomes
        search_handle = Entrez.esearch(db=db, term=term, retmax=retmax, idtype="acc")
        search_results = Entrez.read(search_handle)
        search_handle.close()

        # Get the list of accession numbers
        accession_list = search_results['IdList']

        # Initialize list to store publication info
        publication_info = []

        # Use tqdm to show progress bar
        for acc in tqdm(accession_list, desc="Fetching publications"):
            fetch_handle = Entrez.efetch(db=db, id=acc, rettype="gb", retmode="xml")
            fetch_results = fetch_handle.read()
            fetch_handle.close()
            
            # Parse XML
            root = ET.fromstring(fetch_results)
            publication = root.find(".//GBSeq_references/GBReference")
            
            if publication is not None:
                title = publication.findtext("GBReference_title", default="N/A")
                authors = publication.findtext("GBReference_authors", default="N/A")
                pubdate = publication.findtext("GBReference_journal", default="N/A")
                # Extract year from journal field if possible
                year = pubdate.split()[-1] if pubdate != "N/A" and len(pubdate.split()) > 0 else "N/A"
            else:
                title = "N/A"
                authors = "N/A"
                pubdate = "N/A"
                year = "N/A"
            
            publication_info.append((acc, title, authors, year))
            time.sleep(1)  # Add a short delay between requests to avoid overloading NCBI servers
        
        return publication_info

    except Exception as e:
        print(f"An error occurred: {e}")
        return None

def save_to_csv(filename, data):
    # Write data to CSV file
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        # Write header
        writer.writerow(['Accession', 'Title', 'Authors', 'Publication Year'])
        # Write data rows
        for acc, title, authors, year in data:
            writer.writerow([acc, title, authors, year])

# Example usage
email = "your_email@example.com"
publications = fetch_shigella_genome_accessions_and_publications(email)

if publications is not None:
    # Print the results
    print(f"\nFound {len(publications)} Shigella complete genome sequences with publication details.")
    for acc, title, authors, year in publications:
        print(f"Accession: {acc}, Title: {title}, Authors: {authors}, Publication Year: {year}")

    # Save to CSV file
    csv_filename = "shigella_genome_accessions_with_publications.csv"
    save_to_csv(csv_filename, publications)
    
    print(f"\nResults saved to {csv_filename}")
else:
    print("\nFailed to fetch publications. Check your email and network connection.")

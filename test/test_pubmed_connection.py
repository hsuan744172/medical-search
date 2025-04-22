import sys
import asyncio
from Bio import Entrez
from dotenv import load_dotenv
import os

load_dotenv()
Entrez.email = os.getenv("ENTREZ_EMAIL")
api_key = os.getenv("ENTREZ_API_KEY")
if api_key:
    Entrez.api_key = api_key

print(f"Using email: {Entrez.email}")
print(f"Using API key: {api_key is not None}")

def test_search():
    try:
        print("Testing PubMed search...")
        handle = Entrez.esearch(db="pubmed", term="cancer", retmax=5)
        records = Entrez.read(handle)
        handle.close()
        
        if records['IdList']:
            print(f"Search successful! Found {len(records['IdList'])} results.")
            print(f"IDs: {records['IdList']}")
            
            # Test fetching a paper
            paper_id = records['IdList'][0]
            print(f"\nTesting paper fetch for ID {paper_id}...")
            handle = Entrez.efetch(db="pubmed", id=paper_id, rettype="medline", retmode="xml")
            # Use read instead of parse since the error message specifically mentions this
            paper = Entrez.read(handle)
            handle.close()
            
            if not paper:
                print("No paper data returned.")
                return False
            
            # Print raw structure for debugging
            print(f"Paper structure type: {type(paper)}")
            if hasattr(paper, 'keys'):
                print(f"Available top-level keys: {list(paper.keys())}")
            elif isinstance(paper, list) and paper and hasattr(paper[0], 'keys'):
                print(f"Available keys in first item: {list(paper[0].keys())}")
                
            # Try to access PubmedArticle which is often the container
            if isinstance(paper, dict) and 'PubmedArticle' in paper:
                article = paper['PubmedArticle'][0]
                print(f"Article structure: {type(article)}")
                
                # Try to extract from MedlineCitation -> Article
                if 'MedlineCitation' in article and 'Article' in article['MedlineCitation']:
                    med_article = article['MedlineCitation']['Article']
                    title = med_article.get('ArticleTitle', 'No title')
                    print(f"Title: {title}")
                    
                    # Try to extract year from PubDate
                    year = "Unknown"
                    if 'Journal' in med_article and 'JournalIssue' in med_article['Journal']:
                        if 'PubDate' in med_article['Journal']['JournalIssue']:
                            pub_date = med_article['Journal']['JournalIssue']['PubDate']
                            if 'Year' in pub_date:
                                year = pub_date['Year']
                    print(f"Year: {year}")
                    
                    # Try to extract authors
                    authors = []
                    if 'AuthorList' in med_article and med_article['AuthorList']:
                        for author in med_article['AuthorList']:
                            if 'LastName' in author and 'ForeName' in author:
                                authors.append(f"{author['LastName']} {author['ForeName']}")
                            elif 'LastName' in author:
                                authors.append(author['LastName'])
                            elif 'CollectiveName' in author:
                                authors.append(author['CollectiveName'])
                
                print(f"Authors: {', '.join(authors) if authors else 'Unknown'}")
            else:
                print("Couldn't find expected structure. Raw data:")
                print(paper)
                
            print("\nConnection is working properly!")
        else:
            print("Search returned no results. Check your query.")
    except Exception as e:
        print(f"Error: {str(e)}")
        print(f"Error type: {type(e).__name__}")
        print("\nPossible solutions:")
        print("1. Check your email and API key in .env file")
        print("2. Make sure you have internet connection")
        print("3. Verify your Bio/Entrez installation: pip install biopython")
        return False
    
    return True

if __name__ == "__main__":
    success = test_search()
    if not success:
        sys.exit(1)

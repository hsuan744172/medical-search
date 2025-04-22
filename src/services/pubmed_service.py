from Bio import Entrez
from typing import Dict, List, Any, Union
import logging
import os
from dotenv import load_dotenv
import time

load_dotenv()
Entrez.email = os.getenv("ENTREZ_EMAIL")
if api_key := os.getenv("ENTREZ_API_KEY"):
    Entrez.api_key = api_key

class PubMedService:
    @staticmethod
    def extract_safe_value(data: Any, key: str, default: str = '') -> str:
        """Safely extract values from PubMed data"""
        if key not in data:
            return default
            
        value = data[key]
        if isinstance(value, list) and value:
            return str(value[0])
        elif value is not None:
            return str(value)
        return default

    @staticmethod
    def extract_authors(data: Any) -> List[str]:
        """Safely extract author information"""
        authors = []
        if 'AU' not in data:
            return authors
            
        for author in data['AU']:
            if isinstance(author, dict):
                if 'LastName' in author and 'ForeName' in author:
                    authors.append(f"{author['LastName']} {author['ForeName']}")
                elif 'LastName' in author:
                    authors.append(author['LastName'])
            elif isinstance(author, str):
                authors.append(author)
        return authors

    @staticmethod
    async def search_papers(query: str, max_results: int = 10, retries: int = 3) -> List[Dict]:
        for attempt in range(retries):
            try:
                handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
                records = Entrez.read(handle)
                handle.close()

                if not records['IdList']:
                    return []

                handle = Entrez.efetch(db="pubmed", id=records['IdList'], rettype="medline", retmode="text")
                papers = handle.read().split("\n\n")
                handle.close()
                
                results = []
                for paper_data in papers:
                    if not paper_data.strip():
                        continue
                        
                    lines = paper_data.strip().split('\n')
                    paper = {}
                    current_key = None
                    
                    for line in lines:
                        if not line.strip():
                            continue
                        if line.startswith('PMID-'):
                            paper['PMID'] = line[6:].strip()
                        elif line.startswith('TI  -'):
                            paper['TI'] = [line[6:].strip()]
                            current_key = 'TI'
                        elif line.startswith('AB  -'):
                            paper['AB'] = [line[6:].strip()]
                            current_key = 'AB'
                        elif line.startswith('AU  -'):
                            if 'AU' not in paper:
                                paper['AU'] = []
                            paper['AU'].append(line[6:].strip())
                        elif line.startswith('DP  -'):
                            paper['DP'] = [line[6:].strip()]
                        elif line.startswith('JT  -'):
                            paper['JT'] = [line[6:].strip()]
                        elif line.startswith('TA  -'):
                            paper['TA'] = [line[6:].strip()]
                        # Continue text for multiline fields
                        elif line.startswith('      ') and current_key:
                            paper[current_key][-1] += ' ' + line.strip()
                    
                    paper_dict = {
                        'paperId': paper.get('PMID', ''),
                        'title': paper.get('TI', [''])[0] if paper.get('TI') else '',
                        'abstract': paper.get('AB', [''])[0] if paper.get('AB') else '',
                        'authors': paper.get('AU', []),
                        'year': paper.get('DP', [''])[0][:4] if paper.get('DP') else '',
                        'venue': paper.get('JT', [''])[0] if paper.get('JT') else
                                 paper.get('TA', [''])[0] if paper.get('TA') else '',
                    }
                    results.append(paper_dict)
                
                return results
            except Exception as e:
                logging.error(f"Error searching PubMed (attempt {attempt+1}/{retries}): {str(e)}")
                if attempt < retries - 1:
                    time.sleep(1)  # Wait before retrying
                else:
                    return []

    @staticmethod
    async def fetch_paper_data(paper_id: str, retries: int = 3) -> Dict:
        for attempt in range(retries):
            try:
                handle = Entrez.efetch(db="pubmed", id=paper_id, rettype="medline", retmode="text")
                paper_text = handle.read()
                handle.close()
                
                if not paper_text.strip():
                    raise Exception(f"No paper found with ID {paper_id}")
                
                lines = paper_text.strip().split('\n')
                paper = {}
                current_key = None
                
                for line in lines:
                    if not line.strip():
                        continue
                    if line.startswith('PMID-'):
                        paper['PMID'] = line[6:].strip()
                    elif line.startswith('TI  -'):
                        paper['TI'] = [line[6:].strip()]
                        current_key = 'TI'
                    elif line.startswith('AB  -'):
                        paper['AB'] = [line[6:].strip()]
                        current_key = 'AB'
                    elif line.startswith('AU  -'):
                        if 'AU' not in paper:
                            paper['AU'] = []
                        paper['AU'].append(line[6:].strip())
                    elif line.startswith('DP  -'):
                        paper['DP'] = [line[6:].strip()]
                    elif line.startswith('JT  -'):
                        paper['JT'] = [line[6:].strip()]
                    elif line.startswith('TA  -'):
                        paper['TA'] = [line[6:].strip()]
                    elif line.startswith('MH  -'):
                        if 'MH' not in paper:
                            paper['MH'] = []
                        paper['MH'].append(line[6:].strip())
                    # Continue text for multiline fields
                    elif line.startswith('      ') and current_key:
                        paper[current_key][-1] += ' ' + line.strip()
                
                paper_data = {
                    'title': paper.get('TI', [''])[0] if paper.get('TI') else '',
                    'abstract': paper.get('AB', [''])[0] if paper.get('AB') else '',
                    'authors': paper.get('AU', []),
                    'year': paper.get('DP', [''])[0][:4] if paper.get('DP') else '',
                    'venue': paper.get('JT', [''])[0] if paper.get('JT') else
                            paper.get('TA', [''])[0] if paper.get('TA') else '',
                    'keywords': paper.get('MH', [])
                }
                
                return paper_data
                
            except Exception as e:
                logging.error(f"Error fetching paper data (attempt {attempt+1}/{retries}): {str(e)}")
                if attempt < retries - 1:
                    time.sleep(1)  # Wait before retrying
                else:
                    raise

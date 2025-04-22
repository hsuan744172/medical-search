# Medical Paper Q&A Assistant

A bilingual (English/Chinese) medical research assistant that helps you search and understand medical papers using Google's Gemini AI, PubMed API, and Streamlit.

## Features

- 🔍 **Advanced Paper Search**: Utilize PubMed's API for comprehensive medical paper searches
- 💬 **AI-Powered Q&A**: Get detailed answers about medical papers using Google's Gemini AI
- 🌐 **Bilingual Support**: Seamlessly switch between English and Chinese interfaces
- 📚 **Detailed Paper Analysis**: View comprehensive paper details including abstracts, authors, and metadata
- 💾 **Smart Caching**: Question-answer pairs are cached for improved response times
- 🔄 **Real-time Updates**: Dynamic interface updates with paper selection and search results
- 📊 **Structured Data**: Well-organized paper information with clear categorization

## Project Structure

```
medical-search/
├── src/
│   ├── __init__.py
│   ├── api.py
│   ├── models/
│   │   └── scholar_rag.py
│   └── services/
│       └── pubmed_service.py
├── test/
│   └── test_pubmed_connection.py
├── app.py
├── requirements.txt
├── pyproject.toml
├── .env
└── README.md
```

## Technical Requirements

- Python 3.11+
- Google Cloud Project with:
  - Gemini API enabled
  - Vertex AI TextEmbedding model access (textembedding-gecko@001)
- PubMed API credentials (NCBI Entrez)

## Environment Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/medical-search.git
   cd medical-search
   ```

2. Set up environment:

   **Option 1: Using Poetry (recommended)**
   ```bash
   poetry install
   ```

   **Option 2: Using pip**
   ```bash
   pip install -r requirements.txt
   ```

3. Configure environment variables in .env file:
   ```plaintext
   GOOGLE_API_KEY=your_google_api_key
   ENTREZ_EMAIL=your_email@example.com
   ENTREZ_API_KEY=your_pubmed_api_key
   GOOGLE_CLOUD_PROJECT=your_google_cloud_project_id
   ```

4. Test the PubMed connection:
   ```bash
   python test/test_pubmed_connection.py
   ```

5. Launch the application:
   ```bash
   streamlit run app.py
   ```

## Core Components

- **ScholarRAG**: Core RAG (Retrieval-Augmented Generation) implementation
- **PubMedService**: PubMed API integration service
- **Streamlit Interface**: Bilingual user interface with real-time updates
- **Vector Storage**: FAISS-based vector storage for efficient retrieval
- **Embedding Model**: Google's Vertex AI TextEmbedding model
- **LLM Integration**: Google's Gemini 2.0 Flash model for question answering

## Development Dependencies

Key dependencies:
- streamlit
- google-cloud-aiplatform
- google-generativeai
- langchain & langchain-community
- faiss-cpu
- biopython
- python-dotenv

## Usage Guide

1. **Paper Search**:
   - Use the sidebar search function
   - Enter keywords related to your medical research interest
   - Click "Search" to retrieve relevant papers

2. **Paper Selection**:
   - Browse search results in the sidebar
   - Click on a paper title to select it
   - View detailed paper information

3. **Asking Questions**:
   - With a paper selected, use the chat interface
   - Enter your question about the paper
   - Receive AI-generated responses based on paper content

4. **Language Toggle**:
   - Use the 🌐 button in the header to switch languages
   - Interface dynamically updates to selected language

## Notes

- API rate limits apply for both PubMed (NCBI Entrez) and Google APIs
- Free usage quotas apply to Google's Gemini and TextEmbedding models
- Responses are cached for efficiency
- Paper metadata is stored for quick reference

## License

MIT License
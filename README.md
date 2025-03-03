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
│   ├── main.py
│   ├── models/
│   │   ├── __init__.py
│   │   └── scholar_rag.py
│   └── services/
│       ├── __init__.py
│       └── pubmed_service.py
├── app.py
├── medical_search.py
├── pyproject.toml
└── README.md
```

## Technical Requirements

- Python 3.11+
- Google Cloud Project with:
  - Gemini API enabled
  - Vertex AI TextEmbedding model access(only 90 days free for the model "textembedding-gecko@001")
- PubMed API credentials

## Environment Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/medical-search.git
   cd medical-search
   ```

2. Install dependencies:
   ```bash
   poetry install
   ```

3. Configure environment variables (.env):
   ```plaintext
   GOOGLE_API_KEY=your_google_api_key
   ENTREZ_EMAIL=your_email@example.com
   ENTREZ_API_KEY=your_pubmed_api_key
   ```

4. Set up Google Cloud credentials:
   ```python
   # Update in src/models/scholar_rag.py
   os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = "/path/to/service-account.json"
   ```

5. Launch the application:
   ```bash
   poetry run streamlit run app.py
   ```

## Core Components

- **ScholarRAG**: Core RAG (Retrieval-Augmented Generation) implementation
- **PubMedService**: PubMed API integration service
- **Streamlit Interface**: Bilingual user interface with real-time updates
- **Vector Storage**: FAISS-based vector storage for efficient retrieval
- **Embedding Model**: Google's Vertex AI TextEmbedding model
- **LLM Integration**: Google's Gemini Pro model for question answering

## Development Dependencies

Key dependencies managed through Poetry:
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

- API rate limits apply for both PubMed and Google APIs
- TextEmbedding model has a 90-day free trial period
- Responses are cached for efficiency
- Paper metadata is stored for quick reference

## License

MIT License
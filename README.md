# Full-Context Help Assistant

This application provides an interactive help assistant for cosmology tools like CAMB, GetDist, and Cobaya. It features:

* Support for multiple tools with a tab interface
* Full context documentation for accurate answers
* Multiple AI model options (OpenAI and Gemini)
* Streaming responses for real-time feedback

This application was built using NextJS + ReactJS + TypeScript and supports multiple AI providers.

## Configuration Options

This application supports multiple AI providers and models:

### 1. OpenAI Models

The application uses the OpenAI Chat Completions API with full context documentation. By default, it uses the `gpt-4.1-mini-2025-04-14` model, but you can configure other OpenAI models in the UI.

Benefits:
- Full context documentation for accurate answers
- Streaming responses for real-time feedback
- Model selection in the UI

### 2. Gemini Models

The application also supports Google's Gemini models, including:
- `gemini-2.5-pro-exp-03-25` (Gemini Pro)
- `gemini-2.0-flash` (Gemini Flash)

You can switch between models in the UI based on your needs for speed vs. accuracy.

### Configuration File

The application is configured via the `config.json` file at the root of the project:

```json
{
  "defaultProgram": "camb",
  "showContextLink": true,
  "programs": [
    {
      "id": "camb",
      "name": "CAMB",
      "description": "Code for Anisotropies in the Microwave Background",
      "contextFiles": ["camb-docs.md", "CAMBdemo.md"],
      "combinedContextFile": "camb-combined.md",
      "docsUrl": "https://camb.readthedocs.io/"
    },
    {
      "id": "getdist",
      "name": "GetDist",
      "description": "Plotting and analysis of MCMC samples",
      "contextFiles": ["getdist-readthedocs.md", "plot_gallery.md"],
      "combinedContextFile": "getdist-combined.md",
      "docsUrl": "https://getdist.readthedocs.io/"
    }
  ]
}
```

## Running the Application

1. Ensure you have Node.js installed
2. Clone this repository
3. Install dependencies:
   ```
   npm install
   ```
4. Configure your API keys as environment variables:
   ```
   # For OpenAI
   export OPENAI_API_KEY='your-openai-api-key-here'

   # For Gemini
   export GEMINI_API_KEY='your-gemini-api-key-here'
   ```
5. Build the context files:
   ```
   node scripts/build-context.js && node scripts/generate-context-module.js
   ```
6. Run the development server:
   ```
   npm run dev
   ```
7. Open [http://localhost:3000](http://localhost:3000) in your browser

## Customizing the Application

### Adding New Tools

To add support for additional tools:

1. Update the configuration in `config.json`:
   ```json
   {
     "id": "your-tool-id",
     "name": "Your Tool Name",
     "description": "Description of your tool",
     "contextFiles": ["your-tool-docs.md"],
     "combinedContextFile": "your-tool-combined.md",
     "docsUrl": "https://your-tool-documentation-url/"
   }
   ```

2. Add documentation for your tool in the `context` directory (e.g., `context/your-tool-docs.md`)

3. Run the build scripts to generate the combined context files:
   ```
   node scripts/build-context.js && node scripts/generate-context-module.js
   ```

### Modifying the UI

The UI components are located in the `app/ui` directory. Key components:

- `chat-container.tsx`: Main container that handles program selection and model selection
- `program-tabs.tsx`: Tab interface for switching between different tools
- `chat.tsx`: Implementation of the chat interface
- `model-selector.tsx`: UI for selecting different AI models

### Deployment

This application can be deployed to Vercel with minimal configuration:

1. Connect your GitHub repository to Vercel
2. Configure the environment variables for your API keys
3. Deploy the application

The build process will automatically generate the combined context files during deployment.

## Learn More

This application provides a full-context approach to AI assistance, sending the entire documentation as context to the AI model rather than using RAG (Retrieval Augmented Generation). This approach ensures more accurate and consistent responses, especially for specialized technical domains.

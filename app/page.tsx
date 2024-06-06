import OpenAIAssistant from "@/app/ui/openai-assistant";


export default function Home() {
  return (
    <main>
      <div className="mx-auto mb-12 max-w-lg text-center">
        <div className="m-4">
          <h1 className="mb-4 text-5xl font-extrabold leading-none tracking-tight text-gray-900 md:text-5xl lg:text-5xl">Crossword Compiler Help</h1>
          <div className="mb-6 text-normal font-normal text-gray-500">
          </div>
        </div>
        <OpenAIAssistant 
          assistantId="asst_LozNIVQK0nGDswTpgKlict1B"
          greeting="How can I help you?"
        />
      </div>
    </main>
  );
}

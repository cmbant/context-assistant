'use client'

import { useState, useEffect } from 'react';
import ProgramTabs from './program-tabs';
import Chat from './chat';
import ChatSimple from './chat-simple';
import OpenAIAssistant from './openai-assistant';
import ModelSelector from './model-selector';
import { Program, ChatState, ModelConfig } from '@/app/utils/types';
import { loadConfig, getApiTypeForModel } from '@/app/utils/config';

interface ChatContainerProps {
  programs: Program[];
  defaultProgramId: string;
}

export default function ChatContainer({
  programs,
  defaultProgramId
}: ChatContainerProps) {
  const config = loadConfig();
  const [activeProgram, setActiveProgram] = useState(defaultProgramId);
  const [chatStates, setChatStates] = useState<Record<string, ChatState>>({});
  const [selectedModelId, setSelectedModelId] = useState<string>(config.defaultModelId);

  // Initialize chat states for all programs
  useEffect(() => {
    const initialChatStates: Record<string, ChatState> = {};
    programs.forEach(program => {
      initialChatStates[program.id] = {
        messages: [],
        selectedModelId: selectedModelId
      };
    });
    setChatStates(initialChatStates);
  }, [programs, selectedModelId]);

  // Handle program change
  const handleProgramChange = (programId: string) => {
    setActiveProgram(programId);
  };

  // Handle model change
  const handleModelChange = (modelId: string) => {
    setSelectedModelId(modelId);

    // Update the selected model for all chat states
    setChatStates(prevStates => {
      const newStates = { ...prevStates };
      Object.keys(newStates).forEach(programId => {
        newStates[programId] = {
          ...newStates[programId],
          selectedModelId: modelId
        };
      });
      return newStates;
    });
  };

  // Get the active program
  const getActiveProgram = () => {
    return programs.find(p => p.id === activeProgram) || programs[0];
  };

  const currentProgram = getActiveProgram();
  const showTabs = programs.length > 1;
  const showContextLink = config.showContextLink !== false; // Default to true if not specified
  const apiType = getApiTypeForModel(selectedModelId);

  return (
    <div className="flex flex-col w-full mx-auto border border-gray-300 shadow-lg rounded-md overflow-hidden text-sm sm:text-base">
      <div className="flex flex-col bg-white">
        <div className="flex justify-between items-center">
          <div className="flex-grow">
            {showTabs && (
              <ProgramTabs
                programs={programs}
                activeProgram={activeProgram}
                onProgramChange={handleProgramChange}
              />
            )}
          </div>
          {config.availableModels.length > 1 && (
            <div className="p-2">
              <ModelSelector
                models={config.availableModels}
                selectedModelId={selectedModelId}
                onModelChange={handleModelChange}
              />
            </div>
          )}
        </div>

        <div className="px-4 py-3 text-sm text-gray-700 border-b border-gray-300 bg-white">
          <div className="flex justify-between items-center">
            <div>
              {currentProgram.description}
            </div>
            <div className="flex space-x-4">
              {showContextLink && (
                <a
                  href={`/api/context/${currentProgram.id}`}
                  target="_blank"
                  rel="noopener noreferrer"
                  className="text-blue-500 hover:underline whitespace-nowrap"
                  title="View and download the LLM document context file"
                >
                  <span className="hidden sm:inline">Context</span>
                  <span className="sm:hidden">📄</span>
                </a>
              )}
              <a
                href={currentProgram.docsUrl}
                target="_blank"
                rel="noopener noreferrer"
                className="text-blue-500 hover:underline whitespace-nowrap"
                title="View official documentation"
              >
                <span className="hidden sm:inline">Docs</span>
                <span className="sm:hidden">📖</span>
              </a>
            </div>
          </div>
        </div>
      </div>

      {currentProgram.assistantId ? (
        <OpenAIAssistant
          assistantId={currentProgram.assistantId}
          greeting={`How can I help you?`}
        />
      ) : (
        <ChatSimple
          /* Remove key to prevent re-render when program changes */
          programId={activeProgram}
          greeting={`How can I help you?`}
          apiType={apiType}
          selectedModelId={selectedModelId}
        />
      )}
    </div>
  );
}

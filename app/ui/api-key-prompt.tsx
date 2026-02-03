'use client';

import { useState, useEffect } from 'react';
import { saveUserGeminiKey, getUserGeminiKey, clearUserGeminiKey } from '@/app/utils/usage-tracking';

interface ApiKeyPromptProps {
  isOpen: boolean;
  onClose: () => void;
  onKeySubmit: (key: string) => void;
}

export default function ApiKeyPrompt({ isOpen, onClose, onKeySubmit }: ApiKeyPromptProps) {
  const [apiKey, setApiKey] = useState('');
  const [error, setError] = useState('');
  const [existingKey, setExistingKey] = useState<string | null>(null);

  useEffect(() => {
    if (isOpen) {
      const storedKey = getUserGeminiKey();
      setExistingKey(storedKey);
      if (storedKey) {
        setApiKey(storedKey);
      }
    }
  }, [isOpen]);

  if (!isOpen) return null;

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    const trimmedKey = apiKey.trim();

    if (!trimmedKey) {
      setError('Please enter your API key');
      return;
    }

    // Basic validation - Gemini keys typically start with 'AI'
    if (trimmedKey.length < 20) {
      setError('API key appears to be too short');
      return;
    }

    saveUserGeminiKey(trimmedKey);
    onKeySubmit(trimmedKey);
    setError('');
    onClose();
  };

  const handleClear = () => {
    clearUserGeminiKey();
    setApiKey('');
    setExistingKey(null);
    setError('');
  };

  return (
    <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center z-50">
      <div className="bg-white dark:bg-gray-800 rounded-lg shadow-xl max-w-md w-full mx-4 p-6">
        <h2 className="text-xl font-bold text-gray-900 dark:text-gray-100 mb-4">
          API Key Required
        </h2>

        <p className="text-gray-600 dark:text-gray-300 mb-4">
          You have reached the limit of free API calls. To continue using this
          assistant, please provide your own Gemini API key.
        </p>

        <p className="text-sm text-gray-500 dark:text-gray-400 mb-4">
          You can get a free API key from Google AI Studio:{' '}
          <a
            href="https://ai.google.dev/gemini-api/docs/api-key"
            target="_blank"
            rel="noopener noreferrer"
            className="text-blue-600 dark:text-blue-400 hover:underline"
          >
            Get your API key →
          </a>
        </p>

        <form onSubmit={handleSubmit}>
          <div className="mb-4">
            <label
              htmlFor="api-key"
              className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2"
            >
              Gemini API Key
            </label>
            <input
              type="password"
              id="api-key"
              value={apiKey}
              onChange={(e) => setApiKey(e.target.value)}
              placeholder="Enter your Gemini API key"
              autoComplete="off"
              className="w-full px-3 py-2 border border-gray-300 dark:border-gray-600 rounded-md
                         bg-white dark:bg-gray-700 text-gray-900 dark:text-gray-100
                         focus:outline-none focus:ring-2 focus:ring-blue-500"
            />
            {error && (
              <p className="mt-2 text-sm text-red-600 dark:text-red-400">{error}</p>
            )}
          </div>

          <div className="flex justify-between gap-3">
            <div>
              {existingKey && (
                <button
                  type="button"
                  onClick={handleClear}
                  className="px-4 py-2 text-sm text-red-600 dark:text-red-400
                             hover:text-red-800 dark:hover:text-red-300"
                >
                  Clear Key
                </button>
              )}
            </div>
            <div className="flex gap-3">
              <button
                type="button"
                onClick={onClose}
                className="px-4 py-2 text-sm text-gray-600 dark:text-gray-400
                           hover:text-gray-800 dark:hover:text-gray-200"
              >
                Cancel
              </button>
              <button
                type="submit"
                className="px-4 py-2 bg-blue-600 text-white rounded-md
                           hover:bg-blue-700 focus:outline-none focus:ring-2 focus:ring-blue-500"
              >
                Save Key
              </button>
            </div>
          </div>
        </form>

        <p className="mt-4 text-xs text-gray-500 dark:text-gray-400">
          Your API key is stored locally in your browser and is only used to make calls to Google's Gemini API.
        </p>
      </div>
    </div>
  );
}

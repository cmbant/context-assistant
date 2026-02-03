/**
 * Utility for tracking LLM API calls and managing user API keys
 */

const CALL_COUNT_KEY = 'llm_call_count';
const USER_GEMINI_KEY = 'user_gemini_api_key';

/**
 * Get the current count of LLM calls from localStorage
 * @returns The current call count
 */
export function getCallCount(): number {
  if (typeof window === 'undefined') return 0;
  const count = localStorage.getItem(CALL_COUNT_KEY);
  return count ? parseInt(count, 10) : 0;
}

/**
 * Increment the LLM call count in localStorage
 * @returns The new call count
 */
export function incrementCallCount(): number {
  if (typeof window === 'undefined') return 0;
  const currentCount = getCallCount();
  const newCount = currentCount + 1;
  localStorage.setItem(CALL_COUNT_KEY, newCount.toString());
  return newCount;
}

/**
 * Check if the free LLM calls limit has been exceeded
 * @param maxFreeCalls The maximum number of free calls allowed (0 means no free calls, undefined means unlimited)
 * @returns True if limit is exceeded, false otherwise
 */
export function isFreeLimitExceeded(maxFreeCalls: number | undefined): boolean {
  // undefined means no limit configured - allow unlimited free calls
  if (maxFreeCalls === undefined) return false;
  // 0 means no free calls allowed - immediately require API key
  if (maxFreeCalls === 0) return true;
  // Otherwise check if current count has reached the limit
  return getCallCount() >= maxFreeCalls;
}

/**
 * Get the user's stored Gemini API key
 * @returns The user's Gemini API key or null if not set
 */
export function getUserGeminiKey(): string | null {
  if (typeof window === 'undefined') return null;
  return localStorage.getItem(USER_GEMINI_KEY);
}

/**
 * Save the user's Gemini API key to localStorage
 * @param key The API key to save
 */
export function saveUserGeminiKey(key: string): void {
  if (typeof window === 'undefined') return;
  localStorage.setItem(USER_GEMINI_KEY, key);
}

/**
 * Clear the user's Gemini API key from localStorage
 */
export function clearUserGeminiKey(): void {
  if (typeof window === 'undefined') return;
  localStorage.removeItem(USER_GEMINI_KEY);
}

/**
 * Check if a user API key is required (limit exceeded and no key stored)
 * @param maxFreeCalls The maximum number of free calls allowed
 * @returns True if user needs to provide their own API key
 */
export function requiresUserApiKey(maxFreeCalls: number | undefined): boolean {
  if (!isFreeLimitExceeded(maxFreeCalls)) return false;
  return getUserGeminiKey() === null;
}

/**
 * Check if the user has provided their own API key
 * @returns True if user has their own API key stored
 */
export function hasUserApiKey(): boolean {
  return getUserGeminiKey() !== null;
}

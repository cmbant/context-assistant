'use client'

import { Suspense } from 'react'
import ThemeToggle from './theme-toggle'

// This component is a wrapper for ThemeToggle that uses Suspense
// to prevent hydration issues with useSearchParams in production
export default function ThemeWrapper() {
  return (
    <Suspense fallback={<ThemeFallback />}>
      <ThemeToggle />
    </Suspense>
  )
}

// Simple fallback that looks like the theme toggle button but doesn't do anything
function ThemeFallback() {
  return (
    <div className="p-2 rounded-md bg-gray-100 dark:bg-gray-700 w-9 h-9 flex items-center justify-center">
      <span className="w-5 h-5"></span>
    </div>
  )
}

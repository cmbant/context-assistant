'use client'

import { createContext, useContext, useEffect, useState } from 'react'
import { useSearchParams } from 'next/navigation'

type Theme = 'light' | 'dark'

interface ThemeContextType {
  theme: Theme
  toggleTheme: () => void
}

const ThemeContext = createContext<ThemeContextType>({
  theme: 'light',
  toggleTheme: () => {}
})

export function ThemeProvider({ children }: { children: React.ReactNode }) {
  const [theme, setTheme] = useState<Theme>('light')
  const [mounted, setMounted] = useState(false)
  const searchParams = useSearchParams()

  // Once mounted on client, now we can show the UI
  useEffect(() => {
    setMounted(true)

    // Check for URL parameter first (highest priority)
    const urlTheme = searchParams.get('theme') as Theme | null
    const isValidTheme = urlTheme === 'dark' || urlTheme === 'light'

    // Check for user preference in localStorage
    const savedTheme = localStorage.getItem('theme') as Theme | null

    // Determine which theme to use (URL param > localStorage > system preference)
    let themeToUse: Theme

    if (isValidTheme) {
      // URL parameter takes precedence
      themeToUse = urlTheme as Theme
      // Save to localStorage to persist the choice
      localStorage.setItem('theme', themeToUse)
    } else if (savedTheme) {
      // If no URL param but we have localStorage setting
      themeToUse = savedTheme
    } else {
      // Fall back to system preference
      const systemPrefersDark = window.matchMedia('(prefers-color-scheme: dark)').matches
      themeToUse = systemPrefersDark ? 'dark' : 'light'
    }

    // Apply the selected theme
    setTheme(themeToUse)
    document.documentElement.classList.toggle('dark', themeToUse === 'dark')

    if (themeToUse === 'dark') {
      document.documentElement.style.backgroundColor = '#1a202c'
      document.body.style.backgroundColor = '#1a202c'
    } else {
      document.documentElement.style.backgroundColor = '#ffffff'
      document.body.style.backgroundColor = '#ffffff'
    }
  }, [searchParams])

  const toggleTheme = () => {
    const newTheme = theme === 'light' ? 'dark' : 'light'
    console.log(`Toggling theme from ${theme} to ${newTheme}`)

    // Update state
    setTheme(newTheme)

    // Update classList directly - ensure it's correct
    if (newTheme === 'dark') {
      document.documentElement.classList.add('dark')
      document.documentElement.style.backgroundColor = '#1a202c'
      document.body.style.backgroundColor = '#1a202c'
    } else {
      document.documentElement.classList.remove('dark')
      document.documentElement.style.backgroundColor = '#ffffff'
      document.body.style.backgroundColor = '#ffffff'
    }

    // Save to localStorage
    localStorage.setItem('theme', newTheme)

    console.log('After toggle, classList:', document.documentElement.classList.toString())
  }

  // Avoid hydration mismatch by not rendering anything until mounted
  if (!mounted) {
    return <>{children}</>
  }

  return (
    <ThemeContext.Provider value={{ theme, toggleTheme }}>
      {children}
    </ThemeContext.Provider>
  )
}

export function useTheme() {
  return useContext(ThemeContext)
}

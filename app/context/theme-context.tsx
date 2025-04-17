'use client'

import { createContext, useContext, useEffect, useState } from 'react'

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

  // Once mounted on client, now we can show the UI
  useEffect(() => {
    setMounted(true)

    // Check for user preference in localStorage
    const savedTheme = localStorage.getItem('theme') as Theme | null

    // If user has a saved preference, use that
    if (savedTheme) {
      setTheme(savedTheme)
      document.documentElement.classList.toggle('dark', savedTheme === 'dark')
      if (savedTheme === 'dark') {
        document.documentElement.style.backgroundColor = '#1a202c'
        document.body.style.backgroundColor = '#1a202c'
      } else {
        document.documentElement.style.backgroundColor = '#ffffff'
        document.body.style.backgroundColor = '#ffffff'
      }
    }
    // Otherwise check system preference
    else {
      const systemPrefersDark = window.matchMedia('(prefers-color-scheme: dark)').matches
      setTheme(systemPrefersDark ? 'dark' : 'light')
      document.documentElement.classList.toggle('dark', systemPrefersDark)
      if (systemPrefersDark) {
        document.documentElement.style.backgroundColor = '#1a202c'
        document.body.style.backgroundColor = '#1a202c'
      } else {
        document.documentElement.style.backgroundColor = '#ffffff'
        document.body.style.backgroundColor = '#ffffff'
      }
    }
  }, [])

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

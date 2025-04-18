'use client'

import { createContext, useContext, useEffect, useState } from 'react'

type Theme = 'light' | 'dark'

interface ThemeContextType {
  theme: Theme
  scale: number
  toggleTheme: () => void
}

const ThemeContext = createContext<ThemeContextType>({
  theme: 'light',
  scale: 1,
  toggleTheme: () => {}
})

export function ThemeProvider({ children }: { children: React.ReactNode }) {
  const [theme, setTheme] = useState<Theme>('light')
  const [scale, setScale] = useState<number>(1)
  const [mounted, setMounted] = useState(false)

  // Once mounted on client, now we can show the UI
  useEffect(() => {
    setMounted(true)

    // Function to apply theme
    const applyTheme = (newTheme: Theme) => {
      setTheme(newTheme)
      document.documentElement.classList.toggle('dark', newTheme === 'dark')

      if (newTheme === 'dark') {
        document.documentElement.style.backgroundColor = '#1a202c'
        document.body.style.backgroundColor = '#1a202c'
      } else {
        document.documentElement.style.backgroundColor = '#ffffff'
        document.body.style.backgroundColor = '#ffffff'
      }
    }

    // Function to get URL parameters without using useSearchParams
    const getUrlParam = (name: string): string | null => {
      if (typeof window !== 'undefined') {
        const urlParams = new URLSearchParams(window.location.search)
        return urlParams.get(name)
      }
      return null
    }

    // Check for URL parameters (highest priority)
    const urlTheme = getUrlParam('theme') as Theme | null
    const isValidTheme = urlTheme === 'dark' || urlTheme === 'light'

    // Check for scale parameter
    const urlScale = getUrlParam('scale')
    const scaleValue = urlScale ? parseFloat(urlScale) : null
    const isValidScale = scaleValue !== null && !isNaN(scaleValue) && scaleValue > 0 && scaleValue <= 2

    // Check for user preference in localStorage
    const savedTheme = localStorage.getItem('theme') as Theme | null

    // Check for saved scale in localStorage
    const savedScale = localStorage.getItem('scale')
    const savedScaleValue = savedScale ? parseFloat(savedScale) : null

    // Determine which theme to use (URL param > localStorage > system preference)
    if (isValidTheme) {
      // URL parameter takes precedence
      applyTheme(urlTheme as Theme)
      // Save to localStorage to persist the choice
      localStorage.setItem('theme', urlTheme as Theme)
    } else if (savedTheme) {
      // If no URL param but we have localStorage setting
      applyTheme(savedTheme)
    } else {
      // Fall back to system preference
      const systemPrefersDark = window.matchMedia('(prefers-color-scheme: dark)').matches
      const themeToUse = systemPrefersDark ? 'dark' : 'light'
      applyTheme(themeToUse)
    }

    // Apply scale if valid
    if (isValidScale) {
      setScale(scaleValue as number)
      document.documentElement.style.fontSize = `${scaleValue as number}rem`
      localStorage.setItem('scale', String(scaleValue))
    } else if (savedScaleValue !== null && !isNaN(savedScaleValue) && savedScaleValue > 0 && savedScaleValue <= 2) {
      setScale(savedScaleValue)
      document.documentElement.style.fontSize = `${savedScaleValue}rem`
    } else {
      // Default scale
      document.documentElement.style.fontSize = '1rem'
    }

    // Add an event listener to handle URL changes
    const handleUrlChange = () => {
      // Handle theme changes
      const newUrlTheme = getUrlParam('theme')
      if (newUrlTheme === 'dark' || newUrlTheme === 'light') {
        applyTheme(newUrlTheme as Theme)
        localStorage.setItem('theme', newUrlTheme as Theme)
      }

      // Handle scale changes
      const newUrlScale = getUrlParam('scale')
      const newScaleValue = newUrlScale ? parseFloat(newUrlScale) : null
      if (newScaleValue !== null && !isNaN(newScaleValue) && newScaleValue > 0 && newScaleValue <= 2) {
        setScale(newScaleValue)
        document.documentElement.style.fontSize = `${newScaleValue}rem`
        localStorage.setItem('scale', String(newScaleValue))
      }
    }

    window.addEventListener('popstate', handleUrlChange)
    return () => window.removeEventListener('popstate', handleUrlChange)
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
    <ThemeContext.Provider value={{ theme, scale, toggleTheme }}>
      {children}
    </ThemeContext.Provider>
  )
}

export function useTheme() {
  return useContext(ThemeContext)
}

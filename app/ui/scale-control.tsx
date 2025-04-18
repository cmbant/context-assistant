'use client'

import { useTheme } from '@/app/context/theme-context'
import { useState, useEffect } from 'react'
import { FiZoomIn, FiZoomOut } from 'react-icons/fi'

export default function ScaleControl() {
  const { scale } = useTheme()
  const [currentScale, setCurrentScale] = useState(scale)
  
  // Update local state when context scale changes
  useEffect(() => {
    setCurrentScale(scale)
  }, [scale])
  
  const adjustScale = (newScale: number) => {
    // Ensure scale is within bounds
    newScale = Math.max(0.5, Math.min(2, newScale))
    
    // Update the scale in localStorage
    localStorage.setItem('scale', String(newScale))
    
    // Apply the scale directly
    document.documentElement.style.fontSize = `${newScale}rem`
    
    // Update local state
    setCurrentScale(newScale)
    
    // Update URL parameter without page reload
    const url = new URL(window.location.href)
    url.searchParams.set('scale', String(newScale))
    window.history.replaceState({}, '', url.toString())
  }
  
  const decreaseScale = () => {
    adjustScale(currentScale - 0.1)
  }
  
  const increaseScale = () => {
    adjustScale(currentScale + 0.1)
  }
  
  return (
    <div className="flex items-center space-x-1">
      <button
        onClick={decreaseScale}
        className="p-2 rounded-md bg-gray-100 dark:bg-gray-700 hover:bg-gray-200 dark:hover:bg-gray-600 transition-colors"
        aria-label="Decrease font size"
        title="Decrease font size"
      >
        <FiZoomOut className="w-4 h-4 text-gray-700 dark:text-white" />
      </button>
      
      <span className="text-xs text-gray-500 dark:text-gray-400 w-10 text-center">
        {Math.round(currentScale * 100)}%
      </span>
      
      <button
        onClick={increaseScale}
        className="p-2 rounded-md bg-gray-100 dark:bg-gray-700 hover:bg-gray-200 dark:hover:bg-gray-600 transition-colors"
        aria-label="Increase font size"
        title="Increase font size"
      >
        <FiZoomIn className="w-4 h-4 text-gray-700 dark:text-white" />
      </button>
    </div>
  )
}

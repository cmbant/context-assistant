// This script adds zoom functionality to the page using the viewport meta tag
// It reads the zoom parameter from the URL and applies it to the page

document.addEventListener('DOMContentLoaded', function() {
  try {
    // Function to get URL parameters
    function getUrlParam(name) {
      try {
        if (typeof window !== 'undefined') {
          const urlParams = new URLSearchParams(window.location.search);
          return urlParams.get(name);
        }
      } catch (e) {
        console.error('Error getting URL params:', e);
      }
      return null;
    }

    // Handle zoom parameter
    const urlZoom = getUrlParam('zoom');
    const zoomValue = urlZoom ? parseInt(urlZoom, 10) : null;

    // Only apply valid zoom values (between 50% and 200%)
    if (zoomValue && !isNaN(zoomValue) && zoomValue >= 50 && zoomValue <= 200) {
      // Convert percentage to scale (e.g., 100% -> 1.0, 80% -> 0.8)
      const scale = zoomValue / 100;

      // Get existing viewport meta tag or create a new one
      let viewport = document.querySelector('meta[name="viewport"]');
      if (!viewport) {
        viewport = document.createElement('meta');
        viewport.name = 'viewport';
        document.head.appendChild(viewport);
      }

      // Set the viewport content with the custom scale
      viewport.content = `width=device-width, initial-scale=${scale}`;
      console.log(`Applied zoom scale: ${scale} (${zoomValue}%)`);
    }
  } catch (e) {
    console.error('Error applying zoom:', e);
  }
});

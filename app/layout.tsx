import type { Metadata } from "next";
import { Inter } from "next/font/google";
import "./globals.css";
import "katex/dist/katex.min.css";
import { ThemeProvider } from "./context/theme-context";
import CopyScript from "./ui/copy-script";

const inter = Inter({ subsets: ["latin"] });

export const metadata: Metadata = {
  title: "Cosmology Tools Help Assistant",
  description: "Interactive help assistant for CAMB, GetDist, and other cosmology tools",
};

export default async function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  // Context files should be pre-built and included in the deployment
  // No runtime context generation is needed

  return (
    <html lang="en" suppressHydrationWarning>
      <head>
        <script
          dangerouslySetInnerHTML={{
            __html: `
              (function() {
                // Check for saved theme preference or use browser default
                const theme = localStorage.getItem('theme') ||
                  (window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');

                // Apply theme immediately to prevent flash
                if (theme === 'dark') {
                  document.documentElement.classList.add('dark');
                  document.documentElement.style.backgroundColor = '#1a202c';
                  if (document.body) document.body.style.backgroundColor = '#1a202c';
                } else {
                  document.documentElement.classList.remove('dark');
                  document.documentElement.style.backgroundColor = '#ffffff';
                  if (document.body) document.body.style.backgroundColor = '#ffffff';
                }
              })()
            `,
          }}
        />
        <script
          dangerouslySetInnerHTML={{
            __html: `
              // Add copy buttons to code blocks after the page loads
              document.addEventListener('DOMContentLoaded', function() {
                setTimeout(function() {
                  const preElements = document.querySelectorAll('pre');
                  preElements.forEach(function(pre) {
                    // Make sure the pre element has position relative
                    pre.style.position = 'relative';
                    pre.style.overflow = 'visible';
                  });
                }, 1000);
              });
            `,
          }}
        />
      </head>
      <body className={`${inter.className} min-h-screen`} suppressHydrationWarning>
        <ThemeProvider>
          {children}
          <CopyScript />
        </ThemeProvider>
      </body>
    </html>
  );
}

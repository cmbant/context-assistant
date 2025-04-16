import type { Metadata } from "next";
import { Inter } from "next/font/google";
import "./globals.css";

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
  // In Edge runtime, we can't generate combined context files at startup
  // The combined files should be pre-built and included in the deployment
  // This code is kept for local development compatibility
  try {
    if (typeof window === 'undefined' && process.env.NODE_ENV === 'development') {
      console.log('Running in development mode - context files should be pre-built');
    }
  } catch (error) {
    console.error('Error in layout component:', error);
  }

  return (
    <html lang="en" suppressHydrationWarning>
      <body className={inter.className} suppressHydrationWarning>{children}</body>
    </html>
  );
}

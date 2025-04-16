'use client'

import { useState } from 'react';
import { Program } from '@/app/utils/types';

interface ProgramTabsProps {
  programs: Program[];
  activeProgram: string;
  onProgramChange: (programId: string) => void;
}

export default function ProgramTabs({
  programs,
  activeProgram,
  onProgramChange
}: ProgramTabsProps) {
  return (
    <div className="border-b border-gray-300 bg-gray-100">
      <ul className="flex px-2 pt-2">
        {programs.map((program) => (
          <li key={program.id} className="mr-1">
            <button
              className={`px-5 py-2 rounded-t-md transition-colors ${activeProgram === program.id
                ? 'bg-white text-blue-700 font-semibold border-2 border-gray-300 border-b-0'
                : 'bg-gray-200 text-gray-600 hover:bg-gray-300 border border-gray-300'
              }`}
              onClick={() => onProgramChange(program.id)}
            >
              {program.name}
            </button>
          </li>
        ))}
      </ul>
    </div>
  );
}

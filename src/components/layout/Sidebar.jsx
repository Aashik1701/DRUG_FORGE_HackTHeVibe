import React from 'react';
import { NavLink } from 'react-router-dom';
import {
  LayoutDashboard,
  FlaskConical,
  Layers,
  Dna,
  Target,
  Settings,
  LogOut,
  ChevronLeft,
  ChevronRight,
} from 'lucide-react';
import { useAuth } from '../../context/AuthContext';
import { useDrugForge } from '../../context/DrugForgeContext';

const navItems = [
  { icon: LayoutDashboard, label: 'Dashboard', path: '/app' },
  { icon: FlaskConical, label: 'Lab Bench', path: '/app/analyze' },
  { icon: Target, label: 'Docking Studio', path: '/app/docking' },
  { icon: Dna, label: 'Molecule 3D', path: '/app/visualization' },
  { icon: Layers, label: 'Batch Process', path: '/app/batch' },
  { icon: Settings, label: 'Settings', path: '/app/settings' },
];

const Sidebar = () => {
  const { logout } = useAuth();
  const { isSidebarCollapsed, toggleSidebar } = useDrugForge();

  return (
    <aside
      className={`fixed left-4 top-4 bottom-4 z-50 flex flex-col transition-all duration-300 ease-in-out ${
        isSidebarCollapsed ? 'w-20' : 'w-20 md:w-64'
      }`}
    >
      <div className="relative flex flex-col h-full overflow-hidden border shadow-2xl rounded-3xl border-white/20 bg-white/10 dark:bg-black/30 backdrop-blur-xl">

        {/* Logo Area */}
        <div className="relative flex items-center justify-center h-24 border-b border-white/10">
          <div className="absolute inset-0 opacity-50 bg-gradient-to-r from-bio-teal/10 to-bio-violet/10" />
          {!isSidebarCollapsed ? (
            <h1 className="hidden text-2xl font-bold text-transparent bg-clip-text bg-gradient-to-r from-bio-teal to-bio-violet md:block whitespace-nowrap">
              DrugForge
            </h1>
          ) : null}
          <FlaskConical className={`w-8 h-8 text-bio-teal ${!isSidebarCollapsed ? 'md:hidden' : ''}`} />
        </div>

        {/* Navigation Links */}
        <nav className="flex-1 px-3 py-6 space-y-2 overflow-x-hidden overflow-y-auto">
          {navItems.map((item) => (
            <NavLink
              key={item.path}
              to={item.path}
              end={item.path === '/app'}
              title={isSidebarCollapsed ? item.label : ''}
              className={({ isActive }) => `
                group flex items-center px-3 py-3 rounded-xl transition-all duration-300 relative
                ${isActive
                  ? 'bg-gradient-to-r from-bio-teal/20 to-bio-blue/20 text-bio-teal shadow-lg border border-white/10'
                  : 'text-slate-500 dark:text-slate-400 hover:bg-white/10 hover:text-slate-700 dark:hover:text-slate-200'}
                ${isSidebarCollapsed ? 'justify-center' : ''}
              `}
            >
              {({ isActive }) => (
                <>
                  {isActive && (
                    <div className="absolute top-0 bottom-0 left-0 w-1 rounded-r-full bg-bio-teal" />
                  )}
                  <item.icon className={`w-6 h-6 shrink-0 transition-transform group-hover:scale-110 ${isActive ? 'text-bio-teal' : ''}`} />
                  <span
                    className={`ml-3 font-medium whitespace-nowrap transition-opacity duration-200 ${
                      isSidebarCollapsed ? 'opacity-0 w-0 hidden' : 'opacity-100 hidden md:block'
                    }`}
                  >
                    {item.label}
                  </span>
                </>
              )}
            </NavLink>
          ))}
        </nav>

        {/* Collapse Toggle */}
        <button
          onClick={toggleSidebar}
          className="absolute z-50 items-center justify-center hidden p-1 text-white transition-colors -translate-y-1/2 border rounded-full shadow-lg -right-3 top-1/2 bg-slate-800 border-white/20 hover:bg-bio-teal md:flex"
        >
          {isSidebarCollapsed ? <ChevronRight size={16} /> : <ChevronLeft size={16} />}
        </button>

        {/* Logout */}
        <div className="p-4 border-t border-white/10 bg-black/5">
          <button
            onClick={logout}
            className={`w-full flex items-center ${isSidebarCollapsed ? 'justify-center' : 'justify-center md:justify-start'} px-3 py-3 rounded-xl text-rose-400 hover:bg-rose-500/10 transition-colors`}
          >
            <LogOut className="w-5 h-5 shrink-0" />
            <span className={`ml-3 whitespace-nowrap ${isSidebarCollapsed ? 'hidden' : 'hidden md:block'}`}>Sign Out</span>
          </button>
        </div>
      </div>
    </aside>
  );
};

export default Sidebar;

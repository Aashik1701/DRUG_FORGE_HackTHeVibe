import{j as a}from"./index-CfBLcNG8.js";import"./vendor-D6EHgqcQ.js";import{m as d}from"./utils-BETSRTeM.js";const c=({children:r,className:e="",hoverable:s=!0,hoverEffect:o,animation:n={},onClick:t,...l})=>{const i=o!==void 0?o:s,m={...{initial:{opacity:0,y:20,scale:.95},animate:{opacity:1,y:0,scale:1},exit:{opacity:0,y:-20,scale:.95},transition:{duration:.4,ease:[.23,1,.32,1]}},...n},u=i?"hover:bg-white/20 hover:-translate-y-1 hover:shadow-2xl dark:hover:bg-white/10":"";return a.jsxDEV(d.div,{className:`
        relative overflow-hidden rounded-2xl
        border border-white/20 dark:border-gray-700/30
        bg-white/10 dark:bg-black/20
        backdrop-blur-xl
        shadow-xl
        transition-all duration-300
        ${u}
        ${t?"cursor-pointer":""}
        ${e}
      `,onClick:t,whileHover:i?{y:-5,boxShadow:"0 20px 40px -10px rgba(45, 212, 191, 0.2)"}:{},...m,...l,children:[a.jsxDEV("div",{className:"absolute inset-0 bg-gradient-to-br from-white/5 to-transparent pointer-events-none"},void 0,!1,{fileName:"/Users/mohammedaashik/Documents/PROJECT/AYACOM/Drug_Forge_HackTheVibe/src/components/ui/GlassCard.jsx",lineNumber:57,columnNumber:7},void 0),a.jsxDEV("div",{className:"relative z-10",children:r},void 0,!1,{fileName:"/Users/mohammedaashik/Documents/PROJECT/AYACOM/Drug_Forge_HackTheVibe/src/components/ui/GlassCard.jsx",lineNumber:60,columnNumber:7},void 0),i&&a.jsxDEV("div",{className:"absolute inset-0 opacity-0 hover:opacity-100 transition-opacity duration-300 pointer-events-none",children:a.jsxDEV("div",{className:"absolute inset-0 bg-gradient-to-br from-cyan-500/10 via-transparent to-violet-500/10"},void 0,!1,{fileName:"/Users/mohammedaashik/Documents/PROJECT/AYACOM/Drug_Forge_HackTheVibe/src/components/ui/GlassCard.jsx",lineNumber:67,columnNumber:11},void 0)},void 0,!1,{fileName:"/Users/mohammedaashik/Documents/PROJECT/AYACOM/Drug_Forge_HackTheVibe/src/components/ui/GlassCard.jsx",lineNumber:66,columnNumber:9},void 0)]},void 0,!0,{fileName:"/Users/mohammedaashik/Documents/PROJECT/AYACOM/Drug_Forge_HackTheVibe/src/components/ui/GlassCard.jsx",lineNumber:39,columnNumber:5},void 0)},p=({children:r,className:e="",...s})=>a.jsxDEV(c,{className:`p-8 ${e}`,hoverable:!1,...s,children:r},void 0,!1,{fileName:"/Users/mohammedaashik/Documents/PROJECT/AYACOM/Drug_Forge_HackTheVibe/src/components/ui/GlassCard.jsx",lineNumber:79,columnNumber:5},void 0),f=({children:r,variant:e="primary",className:s="",disabled:o=!1,...n})=>{const t={primary:"bg-cyan-500/20 text-cyan-700 dark:text-cyan-300 border-cyan-500/30 hover:bg-cyan-500/30",secondary:"bg-violet-500/20 text-violet-700 dark:text-violet-300 border-violet-500/30 hover:bg-violet-500/30",success:"bg-emerald-500/20 text-emerald-700 dark:text-emerald-300 border-emerald-500/30 hover:bg-emerald-500/30",danger:"bg-rose-500/20 text-rose-700 dark:text-rose-300 border-rose-500/30 hover:bg-rose-500/30",ghost:"bg-white/10 text-gray-700 dark:text-gray-300 border-white/20 hover:bg-white/20"};return a.jsxDEV(d.button,{className:`
        relative px-6 py-3 rounded-xl
        border backdrop-blur-md
        font-medium transition-all duration-200
        disabled:opacity-50 disabled:cursor-not-allowed
        ${t[e]}
        ${s}
      `,whileHover:o?{}:{scale:1.02,y:-2},whileTap:o?{}:{scale:.98},disabled:o,...n,children:r},void 0,!1,{fileName:"/Users/mohammedaashik/Documents/PROJECT/AYACOM/Drug_Forge_HackTheVibe/src/components/ui/GlassCard.jsx",lineNumber:108,columnNumber:5},void 0)},v=({className:r="",icon:e,...s})=>a.jsxDEV("div",{className:"relative",children:[e&&a.jsxDEV("div",{className:"absolute left-4 top-1/2 -translate-y-1/2 text-gray-500 dark:text-gray-400",children:e},void 0,!1,{fileName:"/Users/mohammedaashik/Documents/PROJECT/AYACOM/Drug_Forge_HackTheVibe/src/components/ui/GlassCard.jsx",lineNumber:134,columnNumber:9},void 0),a.jsxDEV("input",{className:`
          w-full px-4 py-3 rounded-xl
          ${e?"pl-11":""}
          bg-white/30 dark:bg-black/30
          backdrop-blur-md
          border border-white/20 dark:border-gray-700/30
          text-gray-900 dark:text-gray-100
          placeholder-gray-500 dark:placeholder-gray-400
          focus:outline-none focus:ring-2 focus:ring-cyan-500/50 focus:border-cyan-500/50
          transition-all duration-200
          ${r}
        `,...s},void 0,!1,{fileName:"/Users/mohammedaashik/Documents/PROJECT/AYACOM/Drug_Forge_HackTheVibe/src/components/ui/GlassCard.jsx",lineNumber:138,columnNumber:7},void 0)]},void 0,!0,{fileName:"/Users/mohammedaashik/Documents/PROJECT/AYACOM/Drug_Forge_HackTheVibe/src/components/ui/GlassCard.jsx",lineNumber:132,columnNumber:5},void 0),k=({children:r,variant:e="default",className:s=""})=>{const o={default:"bg-white/20 text-gray-700 dark:text-gray-300",primary:"bg-cyan-500/20 text-cyan-700 dark:text-cyan-300",success:"bg-emerald-500/20 text-emerald-700 dark:text-emerald-300",warning:"bg-amber-500/20 text-amber-700 dark:text-amber-300",danger:"bg-rose-500/20 text-rose-700 dark:text-rose-300"};return a.jsxDEV("span",{className:`
      inline-flex items-center px-3 py-1 rounded-full
      text-xs font-medium backdrop-blur-md
      border border-white/20 dark:border-gray-700/30
      ${o[e]}
      ${s}
    `,children:r},void 0,!1,{fileName:"/Users/mohammedaashik/Documents/PROJECT/AYACOM/Drug_Forge_HackTheVibe/src/components/ui/GlassCard.jsx",lineNumber:170,columnNumber:5},void 0)};export{k as G,f as a,c as b,p as c,v as d};

\ProvidesClass{myreport}
\typeout{A class for writing scientific reports based on Article class}
\typeout{Written by Beni Stocker}
\typeout{ }

\LoadClass[10pt,a4paper]{article}

% Language/Charset specific
\RequirePackage[english]{babel}

% Commonly used packages 
\RequirePackage{enumerate}
\RequirePackage{a4wide}
\RequirePackage{amsmath}
\RequirePackage{amssymb}
\RequirePackage{tabularx}
%\RequirePackage{booktabs,url,array}
\RequirePackage[colorlinks=false]{hyperref}
\RequirePackage{pdflscape}
\RequirePackage{caption}
\RequirePackage[pdftex]{graphicx}

% Page layout (margins)
\usepackage[left=3cm,right=5cm]{geometry}

% To use modelsummary() table (dynamically generating model summary tables in R)
\usepackage{tabularray}
\usepackage{float}
\usepackage{graphicx}
\usepackage{codehigh}
\usepackage[normalem]{ulem}
\UseTblrLibrary{booktabs}
\UseTblrLibrary{siunitx}
\newcommand{\tinytableTabularrayUnderline}[1]{\underline{#1}}
\newcommand{\tinytableTabularrayStrikeout}[1]{\sout{#1}}
\NewTableCommand{\tinytableDefineColor}[3]{\definecolor{#1}{#2}{#3}}

% line spacing
%\linespread{1.5}

% % Set running text font
% \RequirePackage[T1]{fontenc}
% \RequirePackage[urw-garamond]{mathdesign}
% \RequirePackage{garamondx}

% Set section heading font
\RequirePackage{titlesec}
\RequirePackage{helvet}
\titleformat{\section}
%  {\normalfont\sffamily\Large\bfseries}
  {\normalfont\sffamily\Large}
  {\thesection}{1em}{}
\titleformat{\subsection}
%  {\normalfont\sffamily\large\bfseries}
  {\normalfont\sffamily\large}
  {\thesubsection}{1em}{}
\titleformat{\subsubsection}
%  {\normalfont\sffamily\large\bfseries}
  {\normalfont\sffamily}
  {\thesubsubsection}{1em}{}

% References
\RequirePackage{natbib}

\RequirePackage{footmisc}

% Nice table horizontal lines
\RequirePackage{booktabs}

\g@addto@macro{\table}{\sffamily}


% Code-formatting
\RequirePackage{verbatim}
\RequirePackage{moreverb}
\let\verbatiminclude\verbatimtabinput
\def\verbatimtabsize{4\relax}

% %%% pdfinclude pdfpages
% \RequirePackage{pdfpages}
% \RequirePackage{pdfpagelabels}
% \RequirePackage{pdftex}
% \RequirePackage{pagebackref}
% \RequirePackage{bookmarks=false}
% \RequirePackage{bookmarksopen}
%\RequirePackage{colorlinks}{hyperref}

%% Redefine \maketitle for Nature style
\def\@maketitle{%
  \newpage\setlength{\parskip}{12pt}%
    {\huge\bfseries\noindent\sloppy \textsf{\@title} \par}%
    %{\noindent\sloppy\textsf{\@date} \par}%
    {\noindent\sloppy\textsf{\@author}}%
}

%% Define the affiliations environment, list each institution as an \item
%% Put after \maketitle
\newenvironment{affiliations}{%
    \setcounter{enumi}{1}%
    \setlength{\parindent}{0in}%
    \slshape\sloppy%
    \begin{list}{\upshape$^{\arabic{enumi}}$}{%
        \usecounter{enumi}%
        \setlength{\leftmargin}{0in}%
        \setlength{\topsep}{0in}%
        \setlength{\labelsep}{0in}%
        \setlength{\labelwidth}{0in}%
        \setlength{\listparindent}{0in}%
        \setlength{\itemsep}{0ex}%
        \setlength{\parsep}{0in}%
        }
    }{\end{list}\par\vspace{12pt}}


\renewenvironment{abstract}{%
\hfill\begin{minipage}{0.95\textwidth}}
{\par\noindent\end{minipage}}

% % math symbols and table lines
% \newcommand{\degree}{\ensuremath{^{\circ}}}
% \newcommand{\tophline}{\hline\noalign{\vspace{1mm}}}
% \newcommand{\middlehline}{\noalign{\vspace{1mm}}\hline\noalign{\vspace{1mm}}}
% \newcommand{\bottomhline}{\noalign{\vspace{1mm}}\hline}
% \newcommand{\PG}[2]{\mbox{\ensuremath{#1\;\mathrm{#2}}}}% Physikalische Groesse

% Captions
%\usepackage[footnotesize]{/home/bstocker/mylatex/sty/caption2}
\renewcommand*{\captionlabelfont}{\small\bfseries}
\renewcommand\captionfont{\small\sffamily}

%\renewcommand{\bibname}{References}

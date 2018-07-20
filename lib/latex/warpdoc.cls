%%
%%
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{warpdoc}
              [2002/10/31 v1.0
 Bernard Parent, 1998-2010]
\newcommand\@ptsize{0}
\newif\if@restonecol
\newif\if@titlepage
\newif\if@Appendix
\newcommand\sectionname{Section}
\@titlepagefalse
\if@compatibility\else
\DeclareOption{a4paper}
   {\setlength\paperheight {297mm}%
    \setlength\paperwidth  {210mm}}
\DeclareOption{a5paper}
   {\setlength\paperheight {210mm}%
    \setlength\paperwidth  {148mm}}
\DeclareOption{b5paper}
   {\setlength\paperheight {250mm}%
    \setlength\paperwidth  {176mm}}
\DeclareOption{letterpaper}
   {\setlength\paperheight {11in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{legalpaper}
   {\setlength\paperheight {14in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{executivepaper}
   {\setlength\paperheight {10.5in}%
    \setlength\paperwidth  {7.25in}}
\DeclareOption{landscape}
   {\setlength\@tempdima   {\paperheight}%
    \setlength\paperheight {\paperwidth}%
    \setlength\paperwidth  {\@tempdima}}
\fi
%\if@compatibility
%  \renewcommand\@ptsize{0}
%\else
%\DeclareOption{10pt}{\renewcommand\@ptsize{0}}
%\fi
%\DeclareOption{11pt}{\renewcommand\@ptsize{1}}
%\DeclareOption{12pt}{\renewcommand\@ptsize{2}}
\@Appendixfalse
\if@compatibility\else
\DeclareOption{oneside}{\@twosidefalse \@mparswitchfalse}
\fi
\DeclareOption{twoside}{\@twosidetrue  \@mparswitchtrue}
\DeclareOption{draft}{\setlength\overfullrule{5pt}}
\if@compatibility\else
\DeclareOption{final}{\setlength\overfullrule{0pt}}
\fi
\DeclareOption{titlepage}{\@titlepagetrue}
\if@compatibility\else
\DeclareOption{notitlepage}{\@titlepagefalse}
\fi
\if@compatibility\else
\DeclareOption{onecolumn}{\@twocolumnfalse}
\fi
\DeclareOption{twocolumn}{\@twocolumntrue}
\DeclareOption{leqno}{\input{leqno.clo}}
\DeclareOption{fleqn}{\input{fleqn.clo}}
\DeclareOption{openbib}{%
  \AtEndOfPackage{%
   \renewcommand\@openbib@code{%
      \advance\leftmargin\bibindent
      \itemindent -\bibindent
      \listparindent \itemindent
      \parsep \z@
      }%
   \renewcommand\newblock{\par}}%
}
\ExecuteOptions{letterpaper,10pt,oneside,onecolumn,final}
\ProcessOptions
\input{size1\@ptsize.clo}
\setlength\lineskip{1\p@}
\setlength\normallineskip{1\p@}
\renewcommand\baselinestretch{}
\setlength\parskip{0\p@ \@plus \p@}
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\setcounter{topnumber}{2}
\renewcommand\topfraction{.7}
\setcounter{bottomnumber}{1}
\renewcommand\bottomfraction{.3}
\setcounter{totalnumber}{3}
\renewcommand\textfraction{.2}
\renewcommand\floatpagefraction{.5}
\setcounter{dbltopnumber}{2}
\renewcommand\dbltopfraction{.7}
\renewcommand\dblfloatpagefraction{.5}


\setcounter{secnumdepth}{3}
\newcounter {part}
\newcounter {section}
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]
\renewcommand\thepart         {\@Roman\c@part}
\renewcommand\thesection      {\@arabic\c@section}
\renewcommand\thesubsection   {\thesection.\@arabic\c@subsection}
\renewcommand\thesubsubsection{\thesubsection .\@arabic\c@subsubsection}
\renewcommand\theparagraph    {\thesubsubsection.\@arabic\c@paragraph}
\renewcommand\thesubparagraph {\theparagraph.\@arabic\c@subparagraph}
\newcommand\part{\par
   \addvspace{4ex}%
   \@afterindentfalse
   \secdef\@part\@spart}

\def\@part[#1]#2{%
    \ifnum \c@secnumdepth >\m@ne
      \refstepcounter{part}%
      \addcontentsline{toc}{part}{\thepart\hspace{1em}#1}%
    \else
      \addcontentsline{toc}{part}{#1}%
    \fi
    {\parindent \z@ \raggedright
     \interlinepenalty \@M
     \normalfont
     \ifnum \c@secnumdepth >\m@ne
       \Large\bfseries\boldmath \partname~\thepart
       \par\nobreak
     \fi
     \huge \bfseries\boldmath #2%
     \markboth{}{}\par}%
    \nobreak
    \vskip 3ex
    \@afterheading}
\def\@spart#1{%
    {\parindent \z@ \raggedright
     \interlinepenalty \@M
     \normalfont
     \huge \bfseries\boltmath #1\par}%
     \nobreak
     \vskip 3ex
     \@afterheading}
\newcommand\section{\@startsection {section}{1}{\z@}%
                                   {-3.5ex \@plus -1ex \@minus -.2ex}%
                                   {2.3ex \@plus.2ex}%
                                   {\normalfont\Large\bfseries\boldmath}}
\newcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {1.5ex \@plus .2ex}%
                                     {\normalfont\large\bfseries\boldmath}}
\newcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {1.5ex \@plus .2ex}%
                                     {\normalfont\normalsize\bfseries\boldmath}}
\newcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                                    {3.25ex \@plus1ex \@minus.2ex}%
                                    {-1em}%
                                    {\normalfont\normalsize\bfseries\boldmath}}
\newcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}%
                                       {3.25ex \@plus1ex \@minus .2ex}%
                                       {-1em}%
                                      {\normalfont\normalsize\bfseries\boldmath}}
\if@twocolumn
  \setlength\leftmargini  {2em}
\else
  \setlength\leftmargini  {2.5em}
\fi
\leftmargin  \leftmargini
\setlength\leftmarginii  {2.2em}
\setlength\leftmarginiii {1.87em}
\setlength\leftmarginiv  {1.7em}
\if@twocolumn
  \setlength\leftmarginv  {.5em}
  \setlength\leftmarginvi {.5em}
\else
  \setlength\leftmarginv  {1em}
  \setlength\leftmarginvi {1em}
\fi
\setlength  \labelsep  {.5em}
\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty
\renewcommand\theenumi{\@arabic\c@enumi}
\renewcommand\theenumii{\@alph\c@enumii}
\renewcommand\theenumiii{\@roman\c@enumiii}
\renewcommand\theenumiv{\@Alph\c@enumiv}
\newcommand\labelenumi{\theenumi.}
\newcommand\labelenumii{(\theenumii)}
\newcommand\labelenumiii{\theenumiii.}
\newcommand\labelenumiv{\theenumiv.}
\renewcommand\p@enumii{\theenumi}
\renewcommand\p@enumiii{\theenumi(\theenumii)}
\renewcommand\p@enumiv{\p@enumiii\theenumiii}
\newcommand\labelitemi{$\m@th\bullet$}
\newcommand\labelitemii{\normalfont\bfseries\boldmath --}
\newcommand\labelitemiii{$\m@th\ast$}
\newcommand\labelitemiv{$\m@th\cdot$}
\newenvironment{description}
               {\list{}{\labelwidth\z@ \itemindent-\leftmargin
                        \let\makelabel\descriptionlabel}}
               {\endlist}
\newcommand*\descriptionlabel[1]{\hspace\labelsep
                                \normalfont\bfseries\boldmath #1}
\newenvironment{verse}
               {\let\\\@centercr
                \list{}{\itemsep      \z@
                        \itemindent   -1.5em%
                        \listparindent\itemindent
                        \rightmargin  \leftmargin
                        \advance\leftmargin 1.5em}%
                \item\relax}
               {\endlist}
\newenvironment{quotation}
               {\list{}{\listparindent 1.5em%
                        \itemindent    \listparindent
                        \rightmargin   \leftmargin
                        \parsep        \z@ \@plus\p@}%
                \item\relax}
               {\endlist}
\newenvironment{quote}
               {\list{}{\rightmargin\leftmargin}%
                \item\relax}
               {\endlist}
\if@compatibility
\newenvironment{titlepage}
    {%
      \if@twocolumn
        \@restonecoltrue\onecolumn
      \else
        \@restonecolfalse\newpage
      \fi
      \thispagestyle{empty}%
      \setcounter{page}\z@
    }%
    {\if@restonecol\twocolumn \else \newpage \fi
    }
\else
\newenvironment{titlepage}
    {%
      \if@twocolumn
        \@restonecoltrue\onecolumn
      \else
        \@restonecolfalse\newpage
      \fi
      \thispagestyle{empty}%
      \setcounter{page}\@ne
    }%
    {\if@restonecol\twocolumn \else \newpage \fi
     \if@twoside\else
        \setcounter{page}\@ne
     \fi
    }
\fi
\newcommand\appendix{\par
  \setcounter{section}{0}%
  \setcounter{subsection}{0}%
  \@Appendixtrue
  \renewcommand\sectionname{Appendix}
  \renewcommand\thesection{\@Alph\c@section}
  \setcounter{equation}{0}
  \@addtoreset{equation}{section}%
  \gdef\theequation{\@Alph\c@section.\arabic{equation}}%
}
\setlength\arraycolsep{5\p@}
\setlength\tabcolsep{6\p@}
\setlength\arrayrulewidth{.4\p@}
\setlength\doublerulesep{2\p@}
\setlength\tabbingsep{\labelsep}
\skip\@mpfootins = \skip\footins
\setlength\fboxsep{3\p@}
\setlength\fboxrule{.4\p@}
\renewcommand\theequation{\@arabic\c@equation}
\newcounter{figure}
\renewcommand\thefigure{\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename~\thefigure}
\newenvironment{figure}
               {\@float{figure}}
               {\end@float}
\newenvironment{figure*}
               {\@dblfloat{figure}}
               {\end@dblfloat}
\newcounter{table}
\renewcommand\thetable{\@arabic\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename~\thetable}
\newenvironment{table}
               {\@float{table}}
               {\end@float}
\newenvironment{table*}
               {\@dblfloat{table}}
               {\end@dblfloat}
\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength\abovecaptionskip{10\p@}
\setlength\belowcaptionskip{0\p@}
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{\fontsizecaption{\sc #1}. #2}%
  \ifdim \wd\@tempboxa >\hsize
    \fontsizecaption{\sc #1}. #2\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries\boldmath}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*\cal{\@fontswitch\relax\mathcal}
\DeclareRobustCommand*\mit{\@fontswitch\relax\mathnormal}
\newcommand\@pnumwidth{1.55em}
\newcommand\@tocrmarg{2.55em}
\newcommand\@dotsep{4.5}
\setcounter{tocdepth}{4}
\newcommand\tableofcontents{%
    \section*{\contentsname
        \@mkboth{%
           \MakeUppercase\contentsname}{\MakeUppercase\contentsname}}%
    \@starttoc{toc}%
    }
\newcommand*\l@part[2]{%
  \ifnum \c@tocdepth >-2\relax
    \addpenalty\@secpenalty
    \addvspace{2.25em \@plus\p@}%
    \begingroup
      \setlength\@tempdima{3em}%
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      {\leavevmode
       \large \bfseries\boldmath #1\hfil \hb@xt@\@pnumwidth{\hss #2}}\par
       \nobreak
       \if@compatibility
         \global\@nobreaktrue
         \everypar{\global\@nobreakfalse\everypar{}}%
      \fi
    \endgroup
  \fi}
\newcommand*\l@section[2]{%
  \ifnum \c@tocdepth >\z@
    \addpenalty\@secpenalty
    \addvspace{1.0em \@plus\p@}%
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode \em
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      #1\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
    \endgroup
  \fi}
\newcommand*\l@subsection{\@dottedtocline{2}{1.5em}{2.3em}}
\newcommand*\l@subsubsection{\@dottedtocline{3}{3.8em}{3.2em}}
\newcommand*\l@paragraph{\@dottedtocline{4}{7.0em}{4.1em}}
\newcommand*\l@subparagraph{\@dottedtocline{5}{10em}{5em}}
\newcommand\listoffigures{%
    \section*{\listfigurename
      \@mkboth{\MakeUppercase\listfigurename}%
              {\MakeUppercase\listfigurename}}%
    \@starttoc{lof}%
    }
\newcommand*\l@figure{\@dottedtocline{1}{1.5em}{2.3em}}
\newcommand\listoftables{%
    \section*{\listtablename
      \@mkboth{%
          \MakeUppercase\listtablename}{\MakeUppercase\listtablename}}%
    \@starttoc{lot}%
    }
\let\l@table\l@figure
\newdimen\bibindent
\setlength\bibindent{1.5em}
\newcommand\newblock{\hskip .11em\@plus.33em\@minus.07em}
\let\@openbib@code\@empty
\newenvironment{theindex}
               {\if@twocolumn
                  \@restonecolfalse
                \else
                  \@restonecoltrue
                \fi
                \columnseprule \z@
                \columnsep 35\p@
                \twocolumn[\section*{\indexname}]%
                \@mkboth{\MakeUppercase\indexname}%
                        {\MakeUppercase\indexname}%
                \thispagestyle{plain}\parindent\z@
                \parskip\z@ \@plus .3\p@\relax
                \let\item\@idxitem}
               {\if@restonecol\onecolumn\else\clearpage\fi}
\newcommand\@idxitem{\par\hangindent 40\p@}
\newcommand\subitem{\@idxitem \hspace*{20\p@}}
\newcommand\subsubitem{\@idxitem \hspace*{30\p@}}
\newcommand\indexspace{\par \vskip 10\p@ \@plus5\p@ \@minus3\p@\relax}
\renewcommand\footnoterule{%
  \kern-3\p@
  \hrule\@width.4\columnwidth
  \kern2.6\p@}
\newcommand\@makefntext[1]{%
    \parindent 1em%
    \noindent
    \hb@xt@1.8em{\hss\@makefnmark}#1}
\newcommand\contentsname{Contents}
\newcommand\listfigurename{List of Figures}
\newcommand\listtablename{List of Tables}
\newcommand\refname{References}
\newcommand\indexname{Index}
\newcommand\figurename{Figure}
\newcommand\tablename{Table}
\newcommand\partname{Part}
\newcommand\appendixname{Appendix}
\newcommand\today{}
\edef\today{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space\number\day, \number\year}
\setlength\columnsep{10\p@}
\setlength\columnseprule{0\p@}
\pagestyle{plain}
\pagenumbering{arabic}
\if@twoside
\else
  \raggedbottom
\fi
\if@twocolumn
  \twocolumn
  \sloppy
  \flushbottom
\else
  \onecolumn
\fi



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fonts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use STIX instead of times: looks almost the same but avoids the bug of times font when used in conjunction with psfrag and a eps fig that also uses times
  \IfFileExists{stix.sty}%
  {
    \AtEndOfClass{
      \RequirePackage[nomath]{stix}
    }
  }
  {
    \AtEndOfClass{
      \RequirePackage{times}  
    }
  }

  \IfFileExists{mtpro2.sty}%
  {
   \AtEndOfClass{
    \RequirePackage{mtpro2}       % mathtimes pro 2 package (commercial)
    \enablesubscriptcorrection
    \DeclareMathSymbol{0}{\mathalpha}{operators}{`0}\DeclareMathSymbol{1}{\mathalpha}{operators}{`1}\DeclareMathSymbol{2}{\mathalpha}{operators}{`2}\DeclareMathSymbol{3}{\mathalpha}{operators}{`3}\DeclareMathSymbol{4}{\mathalpha}{operators}{`4}\DeclareMathSymbol{5}{\mathalpha}{operators}{`5}\DeclareMathSymbol{6}{\mathalpha}{operators}{`6}\DeclareMathSymbol{7}{\mathalpha}{operators}{`7}\DeclareMathSymbol{8}{\mathalpha}{operators}{`8}\DeclareMathSymbol{9}{\mathalpha}{operators}{`9}}
  }
  {
    \AtEndOfClass{
      \RequirePackage{txfonts}
      \RequirePackage{amssymb}
      \let\wtilde\tilde
      \let\wbar\bar
      \let\what\hat
      \let\XL\displaystyle
      \let\mathbold\bf
      \let\overcbrace\overbrace
      \let\undercbrace\underbrace
    }
 }

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Page Layout (Margins and so on)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\RequirePackage{setspace}

%% Footnotes are always single-spaced.
%%\let\@@footnotetext=\@footnotetext
%%\def\@footnotetext{\begin{spacing}{1}\@@footnotetext\end{spacing}}

% Set up page layout.
\setlength{\textheight}{\paperheight} % Height of the main body of the text
\addtolength{\textheight}{-1.4in} 
\setlength{\topmargin}{-0.6in} % 0.7" margin on top of page
%\addtolength{\topmargin}{-\paperheight} 
\setlength{\headsep}{.5in}  % space between header and top of body
\addtolength{\headsep}{-\headheight} % See The LaTeX Companion, p 85
\setlength{\footskip}{.5in}  % space between footer and bottom of body
\setlength{\textwidth}{5.7in} % width of the body of the text
\setlength{\oddsidemargin}{.36in} % 1.25" margin on the left for odd pages
\setlength{\evensidemargin}{0.2in} % 1.25"  margin on the right for even pages

% Marginal notes
\setlength{\marginparwidth}{.75in} % width of marginal notes
\setlength{\marginparsep}{.125in} % space between marginal notes and text

% Make text aligned at the bottom from page to page. comment this out if you
% prefer.
%\flushbottom
% don't force text to be aligned at the bottom
\raggedbottom


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bottom right and bottom center page styles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\ps@botrig{%
   \renewcommand{\@mkboth}{\@gobbletwo}%
   \renewcommand{\@oddhead}{}%
   \renewcommand{\@evenhead}{}%
   \renewcommand{\@evenfoot}{\reset@font\rmfamily\bf\thepage\hfil}%
   \renewcommand{\@oddfoot}{\reset@font\rmfamily\hfil\bf\thepage}}

\def\ps@botcen{%
   \renewcommand{\@mkboth}{\@gobbletwo}%
   \renewcommand{\@oddhead}{}%
   \renewcommand{\@evenhead}{}%
   \renewcommand{\@evenfoot}{\reset@font\rmfamily\hfil\bf\thepage\hfil}%
   \renewcommand{\@oddfoot}{\reset@font\rmfamily\hfil\bf\thepage\hfil}}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fancy header with horizontal line page style
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Define the "headings" page style.
\if@twoside % If two-sided printing.
\def\ps@headings{\let\@mkboth\markboth
   \def\@oddfoot{
      }
   \let\@evenfoot\@oddfoot
   \def\@oddhead{
%      {\begin{picture}(0,0) \put(-3,-5){\line(1,0){377}} \end{picture}} 
       \hspace{-0.2cm} 
       {\small{\sc\rightmark}}
                 \hfil   \rm\it{\thepage} \hspace{-0.15cm}
      }
   \def\@evenhead{
%      {\begin{picture}(0,0) \put(-3,-5){\line(1,0){377}} \end{picture}} 
        \hspace{-0.23cm} 
        \rm\it{\thepage}
                 \hfil  {\small{\sc\leftmark}} \hspace{-0.19cm}
      }
   \def\sectionmark##1{\markboth{\ifnum \c@secnumdepth >\m@ne
      \sectionname ~\thesection. \ \fi ##1}{}}
   \def\subsectionmark##1{\markright{\ifnum \c@secnumdepth >\z@
      \thesubsection. \ \fi ##1}}}
\else % If one-sided printing.
\def\ps@headings{\let\@mkboth\markboth
   \def\@oddfoot{  }
   \def\@oddhead{
%      {\begin{picture}(0,0) \put(-3,-5){\line(1,0){377}} \end{picture}}
       \hspace{-0.2cm} 
            {\small{\sc\rightmark}}
                 \hfil   \rm\it{\thepage} \hspace{-0.15cm}
      }
   \def\sectionmark##1{\markright{\ifnum \c@secnumdepth >\m@ne
      \sectionname ~\thesection. \ \fi ##1}}}
\fi



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Special Groups (TOC, LOF, LOT, LOS, Bibliography)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\RequirePackage{multicol}
\RequirePackage{setspace}

%%%%%%%%%%=-NOMENCLATURE-=%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newlength\nomenclaturelabelwidth
\setlength\nomenclaturelabelwidth{0.13\hsize}
\newlength\nomenclaturecolumnsep
\setlength\nomenclaturecolumnsep{0.09\hsize}

%                                     %
%                                     {1.5ex \@plus .2ex}%
%                                     {\normalfont\normalsize\bfseries\boldmath}}


\newenvironment{nomenclaturelist}[1]
               {\raggedright\vspace{2.25ex}\noindent{\normalfont\normalsize\em #1}\list{}{
	       \itemsep 0 em
	       \topsep  0.6 em
	       \leftmargin \nomenclaturelabelwidth
	       \rightmargin 0.0 in
	       \parsep 0.0 in
	       \itemindent 0.0 in
	       \labelsep 0.06 in
	       \labelwidth \leftmargin
               \addtolength\labelwidth{-\labelsep}
	       \let\makelabel\nomenclaturelistlabel
	       }}
               {\endlist}
\newcommand*\nomenclaturelistlabel[1]{\hspace\labelsep \normalfont #1 \hfill \mbox{=\,}}

\newcommand{\nomenclature}[1]{\gdef\@nomenclature{#1}}


\newcommand{\makenomenclature}
{
 \section*{Nomenclature\@mkboth{Nomenclature}{Nomenclature}}%
% \addcontentsline{toc}{section}{Nomenclature}
 \setlength\columnsep{\nomenclaturecolumnsep}
 \begin{multicols}{2}{
   \small
   \@nomenclature
   \normalsize
 }
 \end{multicols}
}

%%%%%%%%%%=-ABSTRACT-=%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




\newcommand{\abstract}[1]{\gdef\@abstract{#1}}


\newcommand{\makeabstract}
{
 \section*{Abstract\@mkboth{Abstract}{Abstract}}%
 \@abstract
}



%%%%%%%%%%%%%%%%%FANCY title page %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newcommand{\maketitle}
{\begin{titlepage}
   %% Set the line spacing to 1 for the title page.
   \begin{spacing}{1}
   \mbox{}
   \vfill \vfill
   {\begin{picture}(0,0)\thicklines
      %\put(-20,-12){\line(1,0){377}}
      \multiput(-18,-12)(0,0.5){10}{\line(1,0){377}}
    \end{picture}}
   \begin{large}
      \begin{flushright}
         {\Huge \begin{spacing}{1.08} { \bf \@title \\} \end{spacing}}
         ~~\\
         {\huge \emph{\@author}}\\
         ~~\\
         \@date\\
      \end{flushright}
   \end{large}
   {\begin{picture}(0,0) \thicklines
       %\put(-0,45){\line(1,0){377}}
       %\put(0,0){\framebox(377,4){}}
       \multiput(0,35)(0,0.5){5}{\line(1,0){377}}
   \end{picture}}

   \vfill\vfill\vfill\vfill\vfill
   Document generated on \today. \\
   \vfill
   \end{spacing}
   \end{titlepage}
   \newpage
}

\newcommand{\department}[1]{\ifx\empty#1\else\gdef\@department{\ignorespaces#1\unskip}\fi}
\renewcommand{\author}[1]{\ifx\empty#1\else\gdef\@author{\ignorespaces#1\unskip}\fi}
\newcommand{\email}[1]{\ifx\empty#1\else\gdef\@email{\ignorespaces#1\unskip}\fi}
\renewcommand{\title}[1]{\ifx\empty#1\else\gdef\@title{\ignorespaces#1\unskip}\fi}
\newcommand{\institution}[1]{\ifx\empty#1\else\gdef\@institution{\ignorespaces#1\unskip}\fi}

\global\def\@department{Dept.\ of Aerospace Engineering} % Default department is Aerospace
\global\def\@institution{Pusan National University} % default is 
\global\def\@author{\mbox{}}
\global\def\@email{\mbox{}}
\global\def\@title{\mbox{}}

\newcommand{\makewarpdoctitle}{%
 \begin{titlepage}
  \begin{spacing}{1.0}
  \flushright
  ~\\[-1em]
  \vfill
  \begin{minipage}[0in]{\hsize}
    ~~\\
    \center
    \bfseries\boldmath\Huge
    \rule[0.0em]{\hsize}{0.2em}\\[2ex]
    {\@title}\\[2.5ex]
    \Huge
    \rule[0.5em]{\hsize}{0.2em}\\[0.8em]
    \begin{center}
      \large{\@author}\\[0.3em]
      \large{\em\@department,~\@institution}\\[2.2em]
    \end{center}
  \end{minipage}
  \vfill\vfill\vfill\vfill
  \begin{minipage}[0in]{\hsize}
    \parbox{0.72\hsize}{
~\\ ~\\ \vfill\vfill
         \bf \small Written in \@date\\[0.2em]
         \bf \small Copyright \copyright\ \number\year\ by \@author\\[0.2em]
         \bf \small \@email\\[0.2em]
			   \bf \small Document generated on \today.\\ }\hfill
    \parbox{0.23\hsize}{\includegraphics[width=\hsize]{warplogo.eps}}
  \end{minipage}
  \vfill
  \end{spacing}
 \end{titlepage}
}



%%%%%%%%%%=-BIBLIOGRAPHY-=%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newenvironment{thebibliography}[1]
     {\section*{\refname
        \@mkboth{\MakeUppercase\refname}{\MakeUppercase\refname}}%
%      \addcontentsline{toc}{section}{\refname}
      \footnotesize
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \sloppy\clubpenalty4000\widowpenalty4000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}



%%%%%%%%%%%%%%%%FONT SIZES

\newcommand{\fontsizefigure}{\small}
\newcommand{\fontsizetable}{\small}
\newcommand{\fontsizechapterquote}{\small}
\newcommand{\fontsizebibliography}{\small}
\newcommand{\fontsizecaption}{\small}
\newcommand{\fontsizeheader}{\scriptsize}
\newcommand{\fontsizechaptertoc}{\large}


%%%%%%%%%%%%%%%FLOAT PLACEMENT

\setcounter{topnumber}{3}      % maximum number of floats at top of page
\renewcommand\topfraction{.8}  % maximum fraction of page for floats at top
\setcounter{bottomnumber}{3}      % maximum number of floats at page bottom
\renewcommand\bottomfraction{.8}  % (.3) maximum fraction of page for floats at bottom
\setcounter{totalnumber}{5}       % maximum total number of floats on one page
\renewcommand\textfraction{.2}    % minimum fraction of page for text
\renewcommand\floatpagefraction{.8} % minimum fraction of float page that should have floats
\setcounter{dbltopnumber}{3}
\renewcommand\dbltopfraction{.8}
\renewcommand\dblfloatpagefraction{.8}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% \tablecaption command: creates a really nice JCP-like table caption
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newcommand{\tempcaption}{}

\def\tablecaption{\refstepcounter\@captype\@dblarg{\@tablecaption\@captype}}


\long\def\@maketablecaption#1#2{%
    \center\vskip\belowcaptionskip
% put the Table ?? on its own line
      #1 \vspace{0.1em} \par #2\par
    \vskip 0.3em}


\long\def\@tablecaption#1[#2]#3{%
  \baselineskip 14.5\p@
  \addcontentsline{\csname ext@#1\endcsname}{#1}%
    {\protect\numberline{\csname the#1\endcsname}{\ignorespaces #2}}%
  \par
  \begingroup
    \@parboxrestore
    \setbox\@tempboxa\hbox{\csname fnum@#1\endcsname #3}%
      \begingroup
        \let\as@tempcaption\as@shortcaption
        \let\as@shortcaption\as@centercaption
        \@maketablecaption{\fontsizecaption{{ \sc \csname fnum@#1\endcsname.}}}
               {\ignorespaces #3}%
        \let\as@shortcaption\as@tempcaption
      \endgroup
    \par%
  \endgroup}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% table toprule,botrule,midrule from booktab.sty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newdimen\heavyrulewidth
\newdimen\lightrulewidth
\newdimen\cmidrulewidth
\newdimen\belowrulesep
\newdimen\aboverulesep
\newdimen\cmidrulesep
\newdimen\cmidrulekern
\newdimen\defaultaddspace
\heavyrulewidth=.08em
\lightrulewidth=.05em
\cmidrulewidth=.03em
\belowrulesep=.65ex
\aboverulesep=.4ex
\cmidrulesep=\doublerulesep
\cmidrulekern=.25em
\defaultaddspace=.5em
\newcount\rulesflag
\newdimen\@cmidrulewidth
\newcount\@cmidla
\newcount\@cmidlb
\rulesflag=0
\def\toprule{\noalign{\ifnum0=`}\fi
    \@ifnextchar[{\@toprule}{\@toprule[\heavyrulewidth]}}
\def\@toprule[#1]{\ifnum\rulesflag=1\vskip
    \doublerulesep\global\rulesflag=0\fi
    \hrule \@height#1\futurelet\@tempa\@endrule}
\def\@endrule{\ifx\@tempa\toprule\global\rulesflag=1%
    \else\ifx\@tempa\midrule\global\rulesflag=1%
    \else\ifx\@tempa\bottomrule\global\rulesflag=1%
    \else\ifx\@tempa\cmidrule\global\rulesflag=1%
    \else\ifx\@tempa\specialrule\global\rulesflag=1%
    \else\ifx\@tempa\addlinespace\global\rulesflag=1%
    \else\vskip \belowrulesep\fi\fi\fi\fi\fi\fi\ifnum0=`{\fi}}
\def\midrule{\noalign{\ifnum0=`}\fi
    \@ifnextchar[{\@midrule}{\@midrule[\lightrulewidth]}}
\def\@midrule[#1]{\ifnum\rulesflag=1\vskip \doublerulesep\global\rulesflag=0
    \else\vskip \aboverulesep\fi
    \hrule \@height#1\futurelet\@tempa\@endrule}
\def\bottomrule{\noalign{\ifnum0=`}\fi
    \@ifnextchar[{\@midrule}{\@midrule[\heavyrulewidth]}}
\def\addlinespace{\noalign{\ifnum0=`}\fi
    \@ifnextchar[{\@addspace}{\@addspace[\defaultaddspace]}}
\def\@addspace[#1]{\ifnum\rulesflag=1\global\rulesflag=0\fi
    \vskip #1\futurelet\@tempa\@endrule}
\def\specialrule#1#2#3{\noalign{
    \ifnum\rulesflag=1\global\rulesflag=0
    \else\vskip #2\fi\hrule \@height#1\vskip #3}}
\def\cmidrule{\noalign{\ifnum0=`}\fi
    \@ifnextchar[{\@cmidrule}{\@cmidrule[\cmidrulewidth]}}
\def\@cmidrule[#1]{\@ifnextchar({\@@cmidrule[#1]}{\@@cmidrule[#1]()}}
\def\@@cmidrule[#1](#2)#3{\@@@cmidrule[#3]{#1}{#2}}
\def\@@@cmidrule[#1-#2]#3#4{\global\@cmidla#1\relax
    \global\advance\@cmidla\m@ne
    \ifnum\@cmidla>0\global\let\@gtempa\@cmidrulea\else
    \global\let\@gtempa\@cmidruleb\fi
    \global\@cmidlb#2\relax
    \global\advance\@cmidlb-\@cmidla
    \global\@cmidrulewidth=#3
    \global\let\cmlkern@l\z@ \global\let\cmlkern@r\z@
    \@tfor\@tempa :=#4\do{\global\expandafter\let
        \csname cmlkern@\@tempa\endcsname\cmidrulekern}%
\ifnum\rulesflag=0\vskip \aboverulesep\fi\ifnum0=`{\fi}\@gtempa
    \noalign{\ifnum0=`}\fi\futurelet\@tempa\@xcmidrule}
\def\@cmidrulea{\multispan\@cmidla&\multispan\@cmidlb
    \unskip\hskip \cmlkern@l\leaders\hrule \@height\@cmidrulewidth\hfill
    \hskip \cmlkern@r\cr}
\def\@cmidruleb{\multispan\@cmidlb
    \unskip\hskip \cmlkern@l\leaders\hrule \@height\@cmidrulewidth\hfill
    \hskip \cmlkern@r\cr}
\def\@xcmidrule{\ifx\@tempa\cmidrule\vskip-\@cmidrulewidth
    \global\rulesflag=1\else
    \ifx\@tempa\morecmidrules\vskip \cmidrulesep
    \global\rulesflag=1\else
    \vskip \belowrulesep\global\rulesflag=0\fi\fi
    \ifnum0=`{\fi}}
\def\morecmidrules{\noalign{\relax}}

% that's for the spacing in a table activated by the command \addlinespace
\defaultaddspace=.01em


%%%%%%%%%%%%%%%%%%%%%% threeparttable.sty from Donald Arseneau  %%%%%%%%%%%%%%%%%%%%%

% MODIFIED BY BERNARD PARENT, 20 APRIL 2002
% threeparttable.sty   (or 3parttable) (or 3parttab on DOS)
% by Donald Arseneau   Updated on May 13, 1999.
% This file may be distributed, modified, and used in other works with just
% one restriction: modified versions must clearly indicate the modification
% (preferably by a name change).
%
% This package facilitates tables with titles (captions) and notes. The
% title and notes are given a width equal to the body of the table (a
% tabular environment).  By itself, a threeparttable does not float, but
% you can put it in a {table} or a {table*} or some other environment.
% (This causes extra typing, but gives more flexibility.)
%
% Inside a threeparttable there should be a caption, followed by a tabular
% environment, possibly followed by a series of itemized "tablenotes".  At
% present, there is nothing automatic about the notes; you must specify the
% identifier in the body of the table ("\tnote{a}") and in the notes below
% ("\item[a]...").  I chose this method because automatic numbering with
% \footnote would be very hard to use, particularly because many tables
% make repeated reference to a single note.  If someone has a suggestion
% for a convenient, elegant, automatic system, I'll listen!  \tnote
% commands can be given in the caption too, and they will *NOT* appear in
% the list of tables.
%
% \begin{table}
%  \begin{threeparttable}
%   \caption{...}
%   \begin{tabular}...% or {tabular*}
%    ...\tnote{1}&....
%    ...
%   \end{tabular}
%   \begin{tablenotes}
%    \item [1] the first note
%    ...
%   \end{tablenotes}
%  \end{threeparttable}
% \end{table}
%
% NOTE: this version has been altered by BERNARD PARENT to take into
%       account the tablecaption captions
%
\edef\endtablenotes % make a scratch macro for restoring catcodes
  {\catcode\string `\string @=\the\catcode\string`\@
   \catcode\string `\string *=\the\catcode\string`\*}
\catcode`\@=11
\catcode`\*=11

\@ifundefined{@tempboxb}{\@nameuse{newbox}\@tempboxb}{}

\newenvironment{threeparttable}{% 3 parts: title, tabular environment, notes
\relax \ifx\\\@centercr \leavevmode \fi
\vbox\bgroup\sloppy
\@ifundefined{@captype}{\def\@captype{table}}{}%
\topsep\z@
\@enumdepth\z@
\global\let\TPT@hsize\@empty
\let\TPT@LA@caption\@caption \let\@caption\TPT@caption
%%%%%%%% ADDITION BY BERNARD PARENT
\let\TPT@LA@tablecaption\@tablecaption \let\@tablecaption\TPT@tablecaption
%%%%%%%% END OF ADDITION
\TPT@hookin\tabular\TPT@LA@tabular \TPT@hookin\tabular*\TPT@LA@tabular*
\let\TPTrlap\relax}%
{\egroup}

\def\TPT@endtab{\everyhbox{}%
\setbox\@tempboxb\lastbox   % grab tabular environment and measure it
\xdef\TPT@hsize{\hsize\the\wd\@tempboxb \parindent 1em
  \noexpand\@parboxrestore}\TPT@hsize\aftergroup\TPT@hsize
\ifx\TPT@docapt\@undefined\else
  \TPT@docapt \vskip.2\baselineskip
\fi \par \box\@tempboxb \prevdepth\z@ }

\def\TPT@hookin#1#2{\let#2#1\def#1{\edef\TPT@evh{{\the\everyhbox}}\everyhbox
  {\let\tabular\TPT@LA@tabular \let\tabular*\TPT@LA@tabular*
  \aftergroup\TPT@endtab \let\TPTrlap\rlap \everyhbox\TPT@evh}#2}}%

\gdef\TPT@hsize{}

\def\TPT@caption#1[#2]#3{\gdef\TPT@docapt
 {\par\global\let\TPT@docapt\@undefined \TPT@LA@caption{#1}[#2]{#3}}%
 \ifx\TPT@hsize\@empty\else \TPT@docapt \fi \ignorespaces}
%%%%%%%%% ADDITION BY BERNARD PARENT
\def\TPT@tablecaption#1[#2]#3{\gdef\TPT@docapt
 {\par\global\let\TPT@docapt\@undefined \TPT@LA@tablecaption{#1}[#2]{#3}}%
 \ifx\TPT@hsize\@empty\else \TPT@docapt \fi \ignorespaces}
%%%%%%%%% END OF ADDITION
\def\tablenotes{\TPT@hsize \list{}{\topsep\z@skip \partopsep\z@skip
\itemsep\z@ \parsep\z@ \itemindent\z@ \leftmargin1.5em
%  for superscript enumeration, use:
\def\makelabel##1{\tnote{##1}\hfil}\labelwidth .5em\labelsep .3em%
%  for on-the-line enumeration use:
% \labelwidth 1em\labelsep .5em%
\relax}}

\def\tnote#1{\protect\TPTrlap{$^{\fam\z@ #1}$}}%
\def\TPTrlap#1{} % notes in caption will disappear in list of tables!


\endtablenotes % restore catcode of @ and * to starting value

\let\endtablenotes\endlist % and assign proper meaning for \endtablenotes!





\endinput
%%
%% End of file `warpdoc.cls'.

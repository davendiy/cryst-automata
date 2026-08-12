%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         ADM author's package                          %
%                     version 0.2c (June 13, 2003)                      %
%         Copyright (c) 2003 Andrey Popov                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{adm-paper}
             [2003/05/09 v0.2c
                 Algebra and Discrete mathematics paper style]

%\RequirePackage{amsfonts,amsmath,amsthm,amscd,amssymb,latexsym}

\newif\if@abstractsonly

\if@compatibility\else
\DeclareOption{abstractsonly}{\@abstractsonlytrue}
\fi

\DeclareOption{draft}{\PassOptionsToClass{draft}{book}}

\ProcessOptions

\LoadClass[b5paper, 11pt, twoside]{book}


\AtEndDocument{\makepaperend}

%Hooking default section reset
\global\let\section@reset\cl@section

\renewcommand\theequation
   {\@arabic\c@equation}

%Overwriting hyperref commands
\def\texorpdfstring#1#2{#1}

\newcommand{\abstractname}{Abstract}
\renewcommand{\bibname}{References}
\newcommand{\contactname}{Contact information}
\newcommand{\researchname}{RESEARCH ARTICLE}
\newcommand{\surveyname}{SURVEY ARTICLE}
\newcommand{\receivedname}{Received by the editors}

\newcommand\myfootnote[1]{
\renewcommand{\thefootnote}{}
\footnotetext{#1}
\def\thefootnote{\@arabic\c@footnote}
}

\newcommand{\keywordsname}{Key words and phrases}
\newcommand\keywords[1]{
\renewcommand{\thefootnote}{}
\footnotetext{\bf \keywordsname:\space {\it #1.}}
\def\thefootnote{\@arabic\c@footnote}
}
\newcommand{\subjclassname}{Mathematics Subject Classification}
\newcommand\subjclass[2]{
\renewcommand{\thefootnote}{}
\footnotetext{{\bf #1\space\subjclassname:}\space {\it #2.}}
\def\thefootnote{\@arabic\c@footnote}
}

\newlength{\ContentsItemsSep}\setlength{\ContentsItemsSep}{2em}


\setcounter{tocdepth}{0}

\newcounter{papers}
\setcounter{papers}{0}

\def\toclevel@paper{0}

\newcommand{\plaineqnnumbers}
{
\global\let\cl@section\section@reset
\renewcommand\theequation
   {\@arabic\c@equation}
}

\newcommand{\sectioneqnnumbers}
{
\global\let\cl@section\section@reset
\@addtoreset{equation}{section}
\renewcommand\theequation
   {\ifnum \c@section>\z@ \thesection.\fi \@arabic\c@equation}
}

\@addtoreset{figure}{papers}

\renewcommand\thefigure
   {\@arabic\c@figure}

\def\@seccntformat#1{\csname the#1\endcsname.\quad}

\renewcommand \thesection {\@arabic\c@section}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\admcopyright#1{\gdef\@admcopyright{#1}}
\def\@admcopyright{\@latex@warning@no@line{No \noexpand\admcopyright given}}

\def\admtitle#1{\gdef\@admtitle{#1}}
\def\@admtitle{\@latex@warning@no@line{No \noexpand\admtitle given}}

\def\admyear#1{\gdef\@admyear{#1}}
\def\@admyear{\@latex@warning@no@line{No \noexpand\admyear given}}

\def\admnumber#1{\gdef\@admnumber{#1}}
\def\@admnumber{\@latex@warning@no@line{No \noexpand\admnumber given}}

%\def\address#1{\gdef\@address{#1}}
%\def\@address{\@latex@warning@no@line{No \noexpand\addres given}}

\def\title#1{\gdef\@title{#1}}
\def\@title{\@latex@warning@no@line{No \noexpand\title given}}

\def\shorttitle#1{\gdef\@shorttitle{#1}}
\def\@shorttitle{\@latex@warning@no@line{No \noexpand\shorttitle given}}

\def\author#1{\gdef\@author{#1}}
\def\@author{\@latex@warning@no@line{No \noexpand\author given}}

\def\shortauthor#1{\gdef\@shortauthor{#1}}
\def\@shortauthor{\@latex@warning@no@line{No \noexpand\shortauthor given}}

\def\thanks#1{\gdef\@thanks{#1}}
\def\@thanks{\@latex@warning@no@line{No \noexpand\thanks given}}

\def\dedicated#1{\gdef\@dedicated{#1}}
\def\@dedicated{\@latex@warning@no@line{No \noexpand\dedicated given}}


\def\communicated#1{\gdef\@communicated{#1}}
\def\@communicated{\@latex@warning@no@line{No \noexpand\communicated given}}

\def\research{\gdef\@papertype{\researchname}}
\def\survey{\gdef\@papertype{\surveyname}}
\def\@papertype{\@latex@warning@no@line{Paper type is undefined!!!}}


\def\received#1{\gdef\@received{#1}}
\def\@received{\@latex@warning@no@line{No \noexpand\communicated given}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\admcopyright{Journal ``Algebra and Discrete Mathematics''}
\admtitle{Algebra and Discrete Mathematics} \admnumber{?}
%\admvol{I}
\admyear{????}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\global\let\addresses\@empty

\def\address{\@ifnextchar[{\address@optarg}{\address@optarg[\@shortauthor]}}

\def\address@optarg[#1]#2#3#4{
  \typeout{ADDRESS OF: #1}
  \g@addto@macro\addresses{\parbox[t]{10em}{\raggedright \bf #1}\hspace{1em}}
  \g@addto@macro\addresses{\parbox[t]{19em}{#2\par
  %\let\testa\noexpand\csname\string{#4}\endcsname
  %\let\testb\noexpand\csname\string{}\endcsname

  \def\testa{wwwtest#4}
  \let\testb\string{}
  %\@ifundefined{\testa}{\typeout{TEST1}}{\typeout{TEST2}}
  \if\string{\testa}\string{wwwtest}\typeout{OK}\else\typeout{NOT OK}\fi
  \typeout{NEWTEST:\testa,\testb. #4}
  {\it E-Mail:}\space\texttt{#3}\par
  \ifx\testa\testb \typeout{NO URL!!!}\typeout{URL:\testa,\testb.}\typeout{URL: #4}
  \else
  {\it URL:}\space\texttt{#4}\typeout{URL:\testa,\testb. #4}
  \fi
  }\par\vspace{2em}}
}

\global\let\@shorttitle\@empty
\global\let\@shortauthor\@empty
\global\let\@thanks\@empty
\global\let\@dedicated\@empty
\global\let\@communicated\@empty
\global\let\@received\@empty
%\global\let\@papertype\@empty

\newcommand\pagespan{
\if@abstractsonly pp.??--??
\else{\protect\thepage--\pageref{adm\alph{papers}}}
\fi
}

\newcommand\makeadminfo{
{\footnotesize\raggedright
{\@admtitle \ifx\@papertype\@empty \else \hfill \@papertype\par\fi
Number\space\@admnumber.\space{\bfseries(\@admyear).}\space pp.\pagespan\par %P.{\protect\thepage--\pageref{adm\alph{papers}}}\par
\copyright\space\@admcopyright
}
}\par
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MakeTitle definition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand\maketitle{
  %\makepaperend
  \ifx\@shorttitle\@empty \global\let\@shorttitle\@title%\typeout{OK TITLE}
  \fi
  \ifx\@shortauthor\@empty \global\let\@shortauthor\@author%\typeout{OK AUTHOR}
  %\else \typeout{AUTHOR: \@shortauthor}
  \fi
  %\addtocounter{papers}{1}
  \refstepcounter{papers}
  %\let\footnotesize\small
  %\let\footnoterule\relax
  %\let \footnote \thanks
  %\null\vfil
  %\clearpage
  \ifx\@empty\@shorttitle \paper[\@title]{\@title}
  \else \paper[\@shorttitle]{\@title}
  \fi
  \thispagestyle{empty}
  \begin{center}%
    \vskip 1em%
    {\large
     \lineskip .75em%
        {\centering\bfseries \@author}\par}%
      \vskip 1.5em%
    {\ifx\@communicated\@empty\else \small Communicated by\space\@communicated
      \vskip 1.5em\fi}%
    {\ifx\@date\@empty\else\@date \par\fi}%       % Set date in \large size.
    \vskip 0.75em%
    {\it \@dedicated}
  \end{center}\par
  \ifx\@thanks\@empty
  \else\myfootnote{{\it \@thanks}}
  \fi
  \global\let\@thanks\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\@dedicated\@empty
  \global\let\@communicated=\@empty
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand\makepaperend{
\vspace{1cm}
\centerline{\textsc{\contactname}}\medskip
{\raggedright
\addresses}
\ifx\@received\@empty\else{\raggedleft\receivedname:\space\@received.}\fi
\label{adm\alph{papers}}
\global\let\@author\@empty
\global\let\@shortauthor=\@empty
\global\let\@shorttitle=\@empty
\global\let\addresses\@empty
\global\let\@email\@empty
\global\let\@received\@empty
\def\@papertype{\@latex@warning@no@line{Paper type is undefined!!!}}
\global\let\cl@section\section@reset
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



\renewcommand\section{\@startsection {section}{1}{\z@}%
                                   {-3.5ex \@plus -1ex \@minus -.2ex}%
                                   {2.3ex \@plus.2ex}%
                                   {\normalfont\large\bfseries}}
%-------------------------------------------------------------
\newcommand\papername{Paper}

\def\papermark#1{%
      \markboth {#1}{\@shortauthor}}%\typeout{PAPERMARK OK}

\def\sectionmark#1{%
      \ifx\@empty\@shortauthor
      \markright {\@author}%\typeout{MARK SECTION}
      \else
      \markright{\@shortauthor}%\typeout{MARK BY SHORT AUTHOR}
      \fi
      }

\newcommand\@paperapp{\papername}
\newcommand\paper{\clearpage
                    \thispagestyle{plain}%
                    \global\@topnum\z@
                    \@afterindentfalse
                    \secdef\@paper\@spaper}

\newcommand\l@paper[2]{%
  \hbox to\textwidth{#1\hfill{\bf #2}}\vspace{\ContentsItemsSep}}

\newcommand{\contentsfmt}[2]
{\small
\parbox[t]{10em}{\raggedright \bf #1}\hspace{2em}\parbox[t]{19em}{\raggedright #2}%\typeout{TEST CONTENTS}
}

\def\@paper[#1]#2{\ifnum \c@secnumdepth >\m@ne
                       \if@mainmatter
                         \refstepcounter{chapter}%
                         \typeout{\@paperapp\space\thepapers.\space\@author.\space\@title}%
                         %\addcontentsline{toc}{paper}{\protect\contentsfmt{\@shortauthor}{\@title}}\typeout{TEST\@shorttitle}%#1
                         \addcontentsline{toc}{paper}{\texorpdfstring{\protect\contentsfmt{\@shortauthor}{\@title}}{\@author}}
                       \else
                         \addcontentsline{toc}{chapter}{#1}%
                       \fi
                    \else
                      \addcontentsline{toc}{chapter}{#1}%
                    \fi
                    \papermark{#1}%
                    \addtocontents{lof}{\protect\addvspace{10\p@}}%
                    \addtocontents{lot}{\protect\addvspace{10\p@}}%
                    \@makepaperhead{#2}%
                    \@afterheading
                    }


\def\@makepaperhead#1{%
\makeadminfo
  \vspace*{50\p@}%
  {\parindent \z@ \raggedright \normalfont
    \ifnum \c@secnumdepth >\m@ne
    \fi
    \interlinepenalty\@M
    \centering \Large \bfseries #1\par\nobreak
    %\vskip 1\p@
  }}


\def\@spaper#1{
\@makespaperhead{#1}%
\@afterheading
}

\def\@makespaperhead#1{%
  \vspace*{50\p@}%
  {\parindent \z@ \raggedright
    \normalfont
    \interlinepenalty\@M
    \centering \Large \bfseries  #1\par\nobreak
    %\vskip 40\p@
  }}


\newenvironment{abstract}{%
        \small
        \quotation {\scshape\abstractname.}
      }
      {\endquotation
      \if@abstractsonly\makepaperend\endinput\fi
      }

\renewenvironment{thebibliography}[1]{%
  %\@xp\chapter\@xp*\@xp{\bibname}%
\medskip\centerline{{\bf \bibname}}
\normalfont\footnotesize\labelsep .5em\relax
  \renewcommand\theenumiv{\arabic{enumiv}}\let\p@enumiv\@empty
  \list{\@biblabel{\theenumiv}}{\settowidth\labelwidth{\@biblabel{#1}}%
    \leftmargin\labelwidth \advance\leftmargin\labelsep
    \usecounter{enumiv}}%
  \sloppy \clubpenalty\@M \widowpenalty\clubpenalty
  \sfcode`\.=\@m
}{%
  \def\@noitemerr{\@latex@warning{Empty `thebibliography' environment}}%
  \endlist
}

\def\@evenhead{\underline{\vbox{\hbox to \textwidth{\strut\thepage\hfil\scshape\leftmark\hfil}}}}%
\def\@oddhead{\underline{\vbox{\hbox to \textwidth{\hfil\strut{\scshape\rightmark}\hfil\thepage}}}}%


\hoffset=0mm
\voffset=0mm

% top: -25.4mm + 17.5mm = -7.9mm
\topmargin=-7.9mm

%176-35-20
\textwidth=121mm
%250-25-25
\textheight=200mm
% odd: 35-25.4
\oddsidemargin=9.6mm
% even: 20mm-25.4mm
\evensidemargin=-5.4mm

%\newtheorem{theorem}{Theorem}[section]
%\newtheorem{definition}[theorem]{Definition}
%\newtheorem{proposition}[theorem]{Proposition}
%\newtheorem{statement}[theorem]{Statement}
%\newtheorem{lemma}[theorem]{Lemma}
%\newtheorem{corollary}[theorem]{Corollary}
%\newtheorem{example}[theorem]{Example}
%\newtheorem{examples}[theorem]{Examples}
%\newtheorem{remark}[theorem]{Remark}
%\newtheorem{question}[theorem]{Question}

%\newtheorem*{theorem*}{Theorem}
%\newtheorem*{definition*}{Definition}
%\newtheorem*{proposition*}{Proposition}
%\newtheorem*{statement*}{Statement}
%\newtheorem*{lemma*}{Lemma}
%\newtheorem*{corollary*}{Corollary}
%\newtheorem*{example*}{Example}
%\newtheorem*{examples*}{Examples}
%\newtheorem*{remark*}{Remark}
%\newtheorem*{question*}{Question}


\endinput

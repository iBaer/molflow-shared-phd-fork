/*
Program:     MolFlow+
Description: Monte Carlo simulator for ultra-high vacuum
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY
Copyright:   E.S.R.F / CERN
Website:     https://cern.ch/molflow

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
*/
#ifndef _GLFILEBOXH_
#define _GLFILEBOXH_

#include <vector>
#include <time.h>
#include "GLWindow.h"

class GLButton;
class GLList;
class GLLabel;
class GLCombo;

#define MAX_FILTER 256

typedef struct {
  char path[MAX_PATH];
  char file[MAX_PATH];
  char fullName[MAX_PATH];
} FILENAME;

class GLFileBox : public GLWindow {

public:

  static FILENAME *OpenFile(char *path=NULL,char *fileName=NULL,char *title=NULL,const char *filters=NULL,int nbFilter=0);
  static std::vector<FILENAME> OpenMultipleFiles(const char *filters = NULL, const char *title=NULL);
  static FILENAME *SaveFile(char *path=NULL,char *fileName=NULL,char *title=NULL,const char *filters=NULL,int nbFilter=0);

private:

  static void InitDrivePaths();

  GLFileBox(char *path,char *fileName,char *title,char **filters,int nbFilter,bool openMode);
  ~GLFileBox();
  void ProcessMessage(GLComponent *src,int message);
  void UpdateFileList(char *path);
  void Back();
  char *GetSizeStr(DWORD size);
  char *GetTimeStr(time_t time);
  char *GetExtension(char *fileName);
  char *GetName(char *fileName);
  bool  MatchFilters(char *fileName);
  void  AddToPathHist(char *path);
  void  AddToFileHist(char *file);
  bool CheckDirectory(char *dirName);

  GLList      *fileList;
  GLLabel     *filterLabel;
  GLCombo     *filterCombo;
  GLButton    *okButton;
  GLButton    *cancelButton;
  GLButton    *prevButton;
  GLLabel     *pathLabel;
  GLCombo     *pathText;
  GLLabel     *fileLabel;
  GLCombo     *fileText;

  char        *filterName[MAX_FILTER];
  char        *filterExt[MAX_FILTER];
  int          nbFilter;
  int          curFilter;
  char         curPath[MAX_PATH];
  char         curFile[MAX_PATH];
  int          rCode;
  bool         mode;

};

#endif /* _GLFILEBOXH_ */
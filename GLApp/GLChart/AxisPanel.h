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
#ifndef _GLCHARTAXISPANELH_
#define _GLCHARTAXISPANELH_

class GLAxis;
class GLTitledPanel;
class GLButton;
class GLComponent;
class GLCombo;
class GLToggle;
class GLTextField;
class GLTabWindow;
class GLLabel;
class GLChart;

class AxisPanel {

public:

   AxisPanel(GLAxis *a,int axisType,GLChart *parentChart);
   void AddToPanel(GLTabWindow *parent,int pIdx);
   void ProcessMessage(GLComponent *src,int message);

private:

  void commit();
  void error(char *m);

  GLAxis  *pAxis;
  GLChart *pChart;
  int     type;

  GLTitledPanel  *scalePanel;
  GLTitledPanel  *settingPanel;

  GLLabel     *MinLabel;
  GLTextField *MinText;
  GLLabel     *MaxLabel;
  GLTextField *MaxText;
  GLToggle    *AutoScaleCheck;

  GLLabel    *ScaleLabel;
  GLCombo    *ScaleCombo;
  GLToggle   *SubGridCheck;
  GLToggle   *VisibleCheck;
  GLToggle   *OppositeCheck;

  GLCombo    *FormatCombo;
  GLLabel    *FormatLabel;

  GLLabel     *TitleLabel;
  GLTextField *TitleText;

  GLLabel     *ColorLabel;
  GLLabel     *ColorView;
  GLButton    *ColorBtn;

  GLLabel     *PositionLabel;
  GLCombo     *PositionCombo;

};

#endif /* _GLCHARTAXISPANELH_ */

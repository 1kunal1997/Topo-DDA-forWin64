function out = model
%
% DDA_template.m
%
% Model exported on Aug 29 2023, 03:17 by COMSOL 6.1.0.282.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('C:\Users\kl89\Downloads');

model.label('DDA_template.mph');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 3);

model.component('comp1').mesh.create('mesh1');

model.label('DDA_template.mph');

model.component('comp1').geom('geom1').run;

out = model;

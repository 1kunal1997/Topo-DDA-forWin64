function out = model
%
% DDA_geo_cylinder_it300_piecewise_aboslute_binarized_extrude_unionerror.m
%
% Model exported on Feb 20 2024, 22:38 by COMSOL 6.1.0.357.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('E:\Topo-DDA-forWin64\Topo-DDA-forWin64');

model.label('DDA_template.mph');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 3);

model.component('comp1').mesh.create('mesh1');

model.label('DDA_template.mph');

model.component('comp1').geom('geom1').run;

model.label('DDA_template.mph');
model.component.label('Components');
model.label('DDA_geo_cylinder_it300_piecewise_aboslute_binarized_extrude.mph');

model.result.create('pg1', 'PlotGroup3D');
model.result('pg1').label('Electric Field (ewfd)');
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').set('defaultPlotID', 'ElectromagneticWavesFrequencyDomain/phys1/pdef1/pcond1/pg1');
model.result('pg1').feature.create('mslc1', 'Multislice');
model.result('pg1').feature('mslc1').set('smooth', 'internal');
model.result('pg1').feature('mslc1').set('data', 'parent');
model.result.remove('pg1');

out = model;

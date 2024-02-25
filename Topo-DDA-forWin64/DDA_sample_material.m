function out = model
%
% DDA_sample_material.m
%
% Model exported on Nov 14 2023, 20:52 by COMSOL 6.1.0.357.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('E:\Topo-DDA-forWin64\Topo-DDA-forWin64');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 3);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('blk1', 'Block');
model.component('comp1').geom('geom1').feature('blk1').setIndex('layer', 0.5, 0);
model.component('comp1').geom('geom1').runPre('fin');
model.component('comp1').geom('geom1').run;

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').selection.set([2]);
model.component('comp1').material.create('mat2', 'Common');
model.component('comp1').material('mat2').selection.set([1]);

out = model;

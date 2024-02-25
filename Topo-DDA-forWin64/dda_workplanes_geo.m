function out = model
%
% dda_workplanes_geo.m
%
% Model exported on Nov 14 2023, 23:06 by COMSOL 6.1.0.357.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('E:\Topo-DDA-forWin64\Topo-DDA-forWin64');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 3);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('wp1', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp1').set('unite', true);
model.component('comp1').geom('geom1').run('wp1');
model.component('comp1').geom('geom1').create('wp2', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp2').set('unite', true);
model.component('comp1').geom('geom1').feature('wp1').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp2').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').run('wp2');
model.component('comp1').geom('geom1').feature.create('ext1', 'Extrude');
model.component('comp1').geom('geom1').feature('ext1').set('workplane', 'wp1');
model.component('comp1').geom('geom1').feature('ext1').selection('input').set({'wp1'});
model.component('comp1').geom('geom1').run('ext1');
model.component('comp1').geom('geom1').feature.create('ext2', 'Extrude');
model.component('comp1').geom('geom1').feature('ext2').set('workplane', 'wp2');
model.component('comp1').geom('geom1').feature('ext2').selection('input').set({'wp2'});
model.component('comp1').geom('geom1').runPre('fin');

out = model;

function out = model
%
% DDA_geo_fillet.m
%
% Model exported on Nov 14 2023, 17:09 by COMSOL 6.1.0.357.

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
model.label('DDA_geo.mph');

model.component('comp1').geom('geom1').geomRep('cadps');
model.component('comp1').geom('geom1').run('uni1');
model.component('comp1').geom('geom1').create('fil1', 'Fillet3D');
model.component('comp1').geom('geom1').feature('fil1').selection('edge').set('uni1', [437 439]);
model.component('comp1').geom('geom1').feature('fil1').selection('edge').clear('uni1');
model.component('comp1').geom('geom1').feature('fil1').selection('edge').set('uni1', [429 432]);
model.component('comp1').geom('geom1').feature('wp0').geom.run('unil0');
model.component('comp1').geom('geom1').feature('wp0').geom.create('fil1', 'Fillet');
model.component('comp1').geom('geom1').feature('wp0').geom.feature('fil1').selection('pointinsketch').set('unil0', [1 23 34 36 37 39 40 41 42 43 44 45 46 47 48 49 50 54 55 56 69 70 71 75 76 77 78 79 80 81 82 83 84 85 86 88 89 91 102 124]);
model.component('comp1').geom('geom1').feature('wp0').geom.feature('fil1').set('radius', '1[nm]');
model.component('comp1').geom('geom1').feature('wp0').geom.run('fil1');
model.component('comp1').geom('geom1').feature('wp0').geom.feature('fil1').set('radius', '2[nm]');
model.component('comp1').geom('geom1').feature('wp0').geom.run('fil1');
model.component('comp1').geom('geom1').run('ext0');
model.component('comp1').geom('geom1').feature('wp0').geom.feature('fil1').active(false);
model.component('comp1').geom('geom1').feature('wp0').geom.run('fil1');
model.component('comp1').geom('geom1').feature('wp0').geom.run('fil1');
model.component('comp1').geom('geom1').feature('wp0').geom.create('boxsel1', 'BoxSelection');
model.component('comp1').geom('geom1').feature('wp0').geom.feature('boxsel1').set('entitydim', 0);
model.component('comp1').geom('geom1').feature('wp0').geom.feature('boxsel1').set('condition', 'inside');
model.component('comp1').geom('geom1').run('ext0');
model.component('comp1').geom('geom1').create('boxsel1', 'BoxSelection');
model.component('comp1').geom('geom1').feature('wp0').geom.feature.move('boxsel1', 393);
model.component('comp1').geom('geom1').feature('wp0').geom.run('boxsel1');
model.component('comp1').geom('geom1').feature('wp0').geom.selection.create('csel1', 'CumulativeSelection');
model.component('comp1').geom('geom1').feature('wp0').geom.selection('csel1').label('box1');
model.component('comp1').geom('geom1').feature('wp0').geom.feature('boxsel1').set('contributeto', 'csel1');
model.component('comp1').geom('geom1').feature('wp0').geom.feature('fil1').selection('point').named('boxsel1');
model.component('comp1').geom('geom1').feature('wp0').geom.run('fil1');
model.component('comp1').geom('geom1').feature('wp0').geom.feature('fil1').active(true);
model.component('comp1').geom('geom1').feature('wp0').geom.run('fil1');

out = model;

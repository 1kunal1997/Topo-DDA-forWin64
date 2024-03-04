name = "DDA_template.mph";
model = mphload(name);
data = importdata('E:\Topo-DDA-forWin64\Calculations\Random Initial Structure\randomdist_it300_lam542_sym_filterTrue_periodicFalse_beta0_epsilon_0.1_penaltytype_piecewise_absolute0.0to0.5\CoreStructure\CoreStructure299Rounded.txt');

%remove this if you want hybrid material!!!
%bitmask = data;
bitmask = round(data);

Nx = 22;
Ny = 22;
Nz = 10;
dist = 15;

d=strcat(num2str(dist),'[nm]');
h = Nz*dist;
height = strcat(num2str(h),'[nm]');
disp(d);
paraset = configureDictionary("double", "double");
wpnameset = configureDictionary("double", "string");
wpset = dictionary;
%rectsset = configureDictionary("double", "strings");
rectsset = dictionary;
unionlayerset = configureDictionary("double", "string");        % create union of pixels of one type in the work plane to extrude 
unionset = configureDictionary("double", "string");             % create union of material post-extrusion
extrudeset = configureDictionary("double", "string");
selectionset = configureDictionary("double", "string");
selectionlabelset = configureDictionary("double", "string");
materialset = configureDictionary("double", "string");
%disp(set(0));

%model.geom('geom1').run(wpname0);

A = dictionary;
A(0) = {"first"};
A{0}(2:4) = {"second"};
%A(0) = {"first", "second"};
disp(A(0));

% --------------------- Parameters -----------------------------
model.param.set('length', '330[nm]');
model.param.descr('length', '');
model.param.set('width', '330[nm]');
model.param.descr('width', '');
model.param.set('height', '150[nm]');
model.param.descr('height', '');
model.param.set('d', '15[nm]');
model.param.descr('d', 'pixel size');
model.param.set('Vnp', '2.835E-21[m^3]');
model.param.descr('Vnp', '');
model.param.create('par2');
model.param('par2').label('EM Params');
model.param('par2').set('P', '100[uW]');
model.param('par2').descr('P', 'Power of input laser');
model.param('par2').set('w0', '5[um]');
model.param('par2').descr('w0', 'Beam size');
model.param('par2').set('I', 'P/((length/2)*(width/2))');
model.param('par2').descr('I', 'Input intensity');
model.param('par2').set('imp_0', 'sqrt(mu0_const/epsilon0_const)');
model.param('par2').descr('imp_0', 'Vacuum impedence');
model.param('par2').set('n_ext', '1.0');
model.param('par2').descr('n_ext', 'Refractive index of external medium');
model.param('par2').set('eps_ext_re', 'n_ext^2');
model.param('par2').descr('eps_ext_re', 'Real Epsilon of external medium');
model.param('par2').set('eps_ext_im', '0');
model.param('par2').descr('eps_ext_im', 'Imag Epsilon of external medium');
model.param('par2').set('imp_ext', 'sqrt(mu0_const/(epsilon0_const*eps_ext_re))');
model.param('par2').descr('imp_ext', 'Background impedence');
model.param('par2').set('A', 'sqrt(2*I*imp_ext)');
model.param('par2').descr('A', 'Amplitude of input electric field');
model.param('par2').set('lam_min', '500[nm]');
model.param('par2').descr('lam_min', 'Min. lambda');
model.param('par2').set('el_max', 'lam_min/(5*n_ext)');
model.param('par2').descr('el_max', 'Maximum mesh element size');
model.param('par2').set('om_p', '2*pi*8[eV]/h_const');
model.param('par2').descr('om_p', 'Plasma frequency of metal');
model.param('par2').set('gamma', '0.8[eV]/h_const');
model.param('par2').descr('gamma', 'Damping');
model.param('par2').set('lam', '542[nm]');
model.param('par2').descr('lam', 'Lam as param');
model.param('par2').set('f', 'c_const/lam');
model.param('par2').descr('f', 'Freq as param');
model.param('par2').set('eps_inf', '5.5');
model.param('par2').descr('eps_inf', 'High freq permittivity');

% ------------------------- for loop for building workplanes -------------------------------
for j=1:Nx*Ny
    disp(bitmask(3*j));
    para = bitmask(3*j);
    rectname = strcat('rmat', num2str(para*1000), num2str(j));
    rectname = convertCharsToStrings(rectname);
    if isKey(paraset,para)==0           % if this pixel value has not been encountered, create new work plane for it

        paraset(para) = para;
        wpnameset(para) = strcat('workplane', num2str(para*1000));
        model.geom('geom1').create(wpnameset(para), 'WorkPlane');
        wpset(para) = model.geom('geom1').feature(wpnameset(para));
        wpset(para).set('unite', true);
        rectsset(para) = {rectname};
     
    else
        rectsset{para}(length(rectsset{para})+1) = {rectname};

    end
    
    wpset(para).geom().create(rectname, 'Rectangle');
    wpset(para).geom().feature(rectname).set('size', {d,d});
    x = strcat(num2str(mod(j-1,Nx)*dist),'[nm]');
    y = strcat(num2str((ceil(j/Nx)-1)*dist),'[nm]');
   % disp(x);
    wpset(para).geom().feature(rectname).set('pos', {x,y});
    wpset(para).geom().feature(rectname).set('base', 'corner');

end

disp('for loop done!');

% ------------------------- for loop for extrusion and selections -------------------------------
k = values(paraset);
for i = 1:numEntries(paraset)
    key = k(i);
    unionlayerset(key) = strcat('unil', num2str(key*1000));
    wpset(key).geom().create(unionlayerset(key), 'Union');
    wpset(key).geom().feature(unionlayerset(key)).selection('input').set(rectsset{key});
    wpset(key).geom().feature(unionlayerset(key)).set('intbnd', false);
    extrudeset(key) = strcat('extmat', num2str(key*1000));
    model.geom('geom1').create(extrudeset(key), 'Extrude');
    model.geom('geom1').feature(extrudeset(key)).set('workplane', wpnameset(key));
    model.geom('geom1').feature(extrudeset(key)).selection('input').set(wpnameset(key));
    model.geom('geom1').feature(extrudeset(key)).setIndex('distance', height, 0);

    selectionset(key) = strcat('csel', num2str(key*1000));
    selectionlabelset(key) = strcat('Mat', num2str(key*1000));
    model.component('comp1').geom('geom1').selection.create(selectionset(key), 'CumulativeSelection');
    model.component('comp1').geom('geom1').selection(selectionset(key)).label(selectionlabelset(key));
    model.component('comp1').geom('geom1').feature(extrudeset(key)).set('contributeto', selectionset(key));

end

% ------------------- final geometry pieces (cut geom into 1/4, external domain, PML) ---------------------
model.component('comp1').geom('geom1').create('blk1', 'Block');
model.component('comp1').geom('geom1').feature('blk1').set('size', {'width/2' 'length/2' '1'});
model.component('comp1').geom('geom1').feature('blk1').setIndex('size', '10*height', 2);
model.component('comp1').geom('geom1').feature('blk1').set('pos', {'0' '0' '-5*height'});
model.component('comp1').geom('geom1').create('uniall', 'Union');
model.component('comp1').geom('geom1').feature('uniall').selection('input').set(values(extrudeset));
model.component('comp1').geom('geom1').run('uniall');
model.component('comp1').geom('geom1').selection.create('cselall', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('cselall').label('Interior Domain');
model.component('comp1').geom('geom1').feature('uniall').set('contributeto', 'cselall');
model.component('comp1').geom('geom1').create('int1', 'Intersection');
model.component('comp1').geom('geom1').feature('int1').selection('input').set({'blk1' 'uniall'});
model.component('comp1').geom('geom1').feature.duplicate('blk2', 'blk1');
model.component('comp1').geom('geom1').feature('blk2').setIndex('layer', 'height/2', 0);
model.component('comp1').geom('geom1').feature('blk2').set('layertop', true);

model.component('comp1').geom('geom1').runPre('fin');

% -------------- Interpolation Functions for real and imaginary permittivity ---------------------- 
model.component('comp1').func.create('int1', 'Interpolation');
model.component('comp1').func('int1').set('source', 'file');
model.component('comp1').func('int1').set('filename', 'E:\Topo-DDA-forWin64\Topo-DDA-forWin64\diel\TiO2_ALD (raw)Re_eps.txt');
model.component('comp1').func('int1').setIndex('funcs', 'eps_TiO2_re', 0, 0);
model.component('comp1').func('int1').label('Real Permitivitty Titania');
model.component('comp1').func.create('int2', 'Interpolation');
model.component('comp1').func('int2').set('source', 'file');
model.component('comp1').func('int2').set('filename', 'E:\Topo-DDA-forWin64\Topo-DDA-forWin64\diel\TiO2_ALD (raw)Im_eps.txt');
model.component('comp1').func('int2').setIndex('funcs', 'eps_TiO2_im', 0, 0);
model.component('comp1').func('int2').label('Imaginary Permitivitty Titania');

% ----------------- Add physics ---------------------------------
model.component('comp1').physics.create('ewfd', 'ElectromagneticWavesFrequencyDomain', 'geom1');
model.component('comp1').physics('ewfd').feature('wee1').set('DisplacementFieldModel', 'DielectricLoss');

% -------------- Add materials and their properties and assign to selection -------------------
for i = 1:numEntries(paraset)
    key = k(i);
    materialset(key) = strcat('mat', num2str(key*1000));
    model.component('comp1').material.create(materialset(key), 'Common');
    if key == 0
        model.component('comp1').material(materialset(key)).selection.all;
        model.component('comp1').material.move(materialset(key), 0);
    else
        model.component('comp1').material(materialset(key)).selection.named(strcat('geom1_',selectionset(key),'_dom'));
    end
    model.component('comp1').material(materialset(key)).propertyGroup.create('DielectricLoss', 'Dielectric_losses');
    model.component('comp1').material(materialset(key)).propertyGroup('DielectricLoss').set('epsilonBis', {strcat(num2str(key),'*eps_TiO2_im(lam)+',num2str(1-key),'*eps_ext_im')});
    model.component('comp1').material(materialset(key)).propertyGroup('DielectricLoss').set('epsilonPrim', {strcat(num2str(key),'*eps_TiO2_re(lam)+',num2str(1-key),'*eps_ext_re')});
    model.component('comp1').material(materialset(key)).propertyGroup('def').set('relpermeability', {'1'});
    model.component('comp1').material(materialset(key)).propertyGroup('def').set('electricconductivity', {'0'});
   
end

% ----------------------- Create selections for PML, PMC, port -----------------------
model.component('comp1').selection.create('box1', 'Box');
model.component('comp1').selection('box1').set('condition', 'inside');
model.component('comp1').selection('box1').set('zmin', 'height + 1[nm]');
model.component('comp1').selection('box1').label('PML1');
model.component('comp1').selection.duplicate('box2', 'box1');
model.component('comp1').selection('box2').label('PML2');
model.component('comp1').selection('box2').set('zmin', -Inf);
model.component('comp1').selection('box2').set('zmax', '-1[nm]');
model.component('comp1').selection.create('uni1', 'Union');
model.component('comp1').selection('uni1').label('PMLs');
model.component('comp1').selection('uni1').set('input', {'box1' 'box2'});
model.component('comp1').selection.create('box3', 'Box');
model.component('comp1').selection('box3').set('entitydim', 2);
model.component('comp1').selection('box3').set('condition', 'inside');
model.component('comp1').selection('box3').label('PMC1');
model.component('comp1').selection('box3').set('ymin', '-0.5*d');
model.component('comp1').selection('box3').set('ymax', '0.5*d');
model.component('comp1').selection.duplicate('box4', 'box3');
model.component('comp1').selection('box4').label('PMC2');
model.component('comp1').selection('box4').set('ymin', '-0.5*d + width/2');
model.component('comp1').selection('box4').set('ymax', '0.5*d + width/2');
model.component('comp1').selection.create('uni2', 'Union');
model.component('comp1').selection('uni2').set('entitydim', 2);
model.component('comp1').selection('uni2').set('input', {'box3' 'box4'});
model.component('comp1').selection('uni2').label('PMCs');

model.component('comp1').selection.create('box5', 'Box');
model.component('comp1').selection('box5').label('Ext1');
model.component('comp1').selection('box5').set('zmin', '-1[nm]');
model.component('comp1').selection('box5').set('zmax', '-1[nm]');
model.component('comp1').selection('box5').set('condition', 'intersects');
model.component('comp1').selection.duplicate('box6', 'box5');
model.component('comp1').selection('box6').label('Ext2');
model.component('comp1').selection('box6').set('zmin', 'height+1[nm]');
model.component('comp1').selection('box6').set('zmax', 'height+1[nm]');
model.component('comp1').selection.create('uni3', 'Union');
model.component('comp1').selection('uni3').set('input', {'box5' 'box6'});
model.component('comp1').selection('uni3').label('External Domain');

model.component('comp1').selection.create('box7', 'Box');
model.component('comp1').selection('box7').label('Port');
model.component('comp1').selection('box7').set('entitydim', 2);
model.component('comp1').selection('box7').set('condition', 'inside');
model.component('comp1').selection('box7').set('zmin', '5*height - height/2');
model.component('comp1').selection('box7').set('zmax', '5*height - height/2 + 1[nm]');

model.component('comp1').selection.create('box8', 'Box');
model.component('comp1').selection('box8').set('entitydim', 2);
model.component('comp1').selection('box8').set('condition', 'inside');
model.component('comp1').selection('box8').set('zmin', 0);
model.component('comp1').selection('box8').set('zmax', '1[nm]');
model.component('comp1').selection('box8').label('Workplane');

% ------------------------- Variables --------------------------------
model.component('comp1').cpl.create('aveop1', 'Average');
model.component('comp1').cpl('aveop1').set('axisym', true);
model.component('comp1').cpl('aveop1').selection.named('geom1_csel1000_dom');
model.component('comp1').cpl('aveop1').label('Average in Pure Titania');
model.component('comp1').cpl.create('intop1', 'Integration');
model.component('comp1').cpl('intop1').set('axisym', true);
model.component('comp1').cpl('intop1').selection.named('geom1_csel1000_dom');
model.component('comp1').cpl('intop1').label('Integral over Pure Titania');

model.component('comp1').variable.create('var1');
model.component('comp1').variable('var1').set('avE', 'aveop1(ewfd.normE)/A');

model.component('comp1').coordSystem.create('pml1', 'PML');
model.component('comp1').coordSystem('pml1').selection.named('uni1');

% ------------------------ Physics setup -----------------------------
model.component('comp1').physics('ewfd').create('pmc1', 'PerfectMagneticConductor', 2);
model.component('comp1').physics('ewfd').feature('pmc1').selection.named('uni2');
model.component('comp1').physics('ewfd').prop('EquationForm').setIndex('form', 'Frequency', 0);
model.component('comp1').physics('ewfd').prop('EquationForm').setIndex('freq_src', 'userdef', 0);
model.component('comp1').physics('ewfd').prop('EquationForm').setIndex('freq', 'f', 0);
model.component('comp1').physics('ewfd').create('port1', 'Port', 2);
model.component('comp1').physics('ewfd').feature('port1').selection.named('box7');
model.component('comp1').physics('ewfd').feature('port1').set('Pin', 'P');
model.component('comp1').physics('ewfd').feature('port1').set('PortSlit', true);
model.component('comp1').physics('ewfd').feature('port1').set('SlitType', 'DomainBacked');
model.component('comp1').physics('ewfd').feature('port1').set('E0', [1 0 0]);
model.component('comp1').physics('ewfd').feature('port1').set('beta', 'ewfd.k0');

% --------------------- Mesh ------------------------------------
model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri1').selection.named('box8');
model.component('comp1').mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('custom', true);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hmax', 'd/3');
model.component('comp1').mesh('mesh1').create('swe1', 'Sweep');
model.component('comp1').mesh('mesh1').feature('swe1').selection.geom('geom1', 3);
model.component('comp1').mesh('mesh1').feature('swe1').selection.named('uni3');
model.component('comp1').mesh('mesh1').feature('swe1').selection.named('geom1_cselall_dom');
model.component('comp1').mesh('mesh1').feature('swe1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('swe1').feature('size1').set('custom', true);
model.component('comp1').mesh('mesh1').feature('swe1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('swe1').feature('size1').set('hmax', '2*d');
model.component('comp1').mesh('mesh1').create('ftet1', 'FreeTet');
model.component('comp1').mesh('mesh1').feature('ftet1').selection.geom('geom1', 3);
model.component('comp1').mesh('mesh1').feature('ftet1').selection.named('uni3');
model.component('comp1').mesh('mesh1').feature('ftet1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('custom', true);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hmax', 'el_max');
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hmin', 'd/3');
model.component('comp1').mesh('mesh1').create('swe2', 'Sweep');
model.component('comp1').mesh('mesh1').feature('swe2').create('dis1', 'Distribution');
model.component('comp1').mesh('mesh1').run('swe2');

% ----------------------------- Create Study --------------------------------
model.study.create('std1');
model.study('std1').create('stat', 'Stationary');
model.study('std1').feature('stat').set('useparam', true);
model.study('std1').feature('stat').setIndex('pname', 'A', 0);
model.study('std1').feature('stat').setIndex('plistarr', '', 0);
model.study('std1').feature('stat').setIndex('punit', 'V/m', 0);
model.study('std1').feature('stat').setIndex('pname', 'A', 0);
model.study('std1').feature('stat').setIndex('plistarr', '', 0);
model.study('std1').feature('stat').setIndex('punit', 'V/m', 0);
model.study('std1').feature('stat').setIndex('pname', 'lam', 0);
model.study('std1').feature('stat').setIndex('plistarr', 'range(lam_min,1[nm],560[nm])', 0);
model.study('std1').feature('stat').setIndex('punit', 'nm', 0);

model.sol.create('sol1');
model.sol('sol1').study('std1');

model.study('std1').feature('stat').set('notlistsolnum', 1);
model.study('std1').feature('stat').set('notsolnum', 'auto');
model.study('std1').feature('stat').set('listsolnum', 1);
model.study('std1').feature('stat').set('solnum', 'auto');

model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'stat');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'stat');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').set('stol', 0.001);
model.sol('sol1').feature('s1').create('p1', 'Parametric');
model.sol('sol1').feature('s1').feature.remove('pDef');
model.sol('sol1').feature('s1').feature('p1').set('control', 'stat');
model.sol('sol1').feature('s1').set('control', 'stat');
model.sol('sol1').feature('s1').feature('aDef').set('complexfun', true);
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').create('i1', 'Iterative');
model.sol('sol1').feature('s1').feature('i1').set('linsolver', 'bicgstab');
model.sol('sol1').feature('s1').feature('i1').set('prefuntype', 'left');
model.sol('sol1').feature('s1').feature('i1').label('Suggested Iterative Solver (ewfd)');
model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').create('sv1', 'SORVector');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sv1').set('prefun', 'sorvec');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sv1').set('iter', 2);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sv1').set('relax', 1);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sv1').set('sorvecdof', {'comp1_E'});
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').create('sv1', 'SORVector');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sv1').set('prefun', 'soruvec');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sv1').set('iter', 2);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sv1').set('relax', 1);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sv1').set('sorvecdof', {'comp1_E'});
model.sol('sol1').feature('s1').create('i2', 'Iterative');
model.sol('sol1').feature('s1').feature('i2').set('linsolver', 'fgmres');
model.sol('sol1').feature('s1').feature('i2').label('Suggested Iterative Solver (ewfd) 2');
model.sol('sol1').feature('s1').feature('i2').create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').create('sv1', 'SORVector');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').feature('sv1').set('prefun', 'sorvec');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').feature('sv1').set('iter', 2);
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').feature('sv1').set('relax', 1);
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').feature('sv1').set('sorvecdof', {'comp1_E'});
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').create('sv1', 'SORVector');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').feature('sv1').set('prefun', 'soruvec');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').feature('sv1').set('iter', 2);
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').feature('sv1').set('relax', 1);
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').feature('sv1').set('sorvecdof', {'comp1_E'});
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'i1');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').attach('std1');
model.sol('sol1').feature('s1').feature('dDef').active(true);


model.geom.run();
disp("im done");

model.save('E:\Topo-DDA-forWin64\Calculations\Random Initial Structure\randomdist_it300_lam542_sym_filterTrue_periodicFalse_beta0_epsilon_0.1_penaltytype_piecewise_absolute0.0to0.5\binarizedit299.mph')
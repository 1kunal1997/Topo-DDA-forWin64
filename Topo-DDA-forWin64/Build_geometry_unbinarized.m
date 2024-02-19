name = "DDA_template.mph";
model = mphload(name);
data = importdata('cylinder_it300_lam542_sym_epsilon_0.1_penaltytype_piecewise_absolute_0.0to0.5_repeat\CoreStructure\CoreStructure299Rounded.txt');
bitmask = data;
Nx = 22;
Ny = 22;
Nz = 10;
dist = 15;
d=strcat(num2str(dist),'[nm]');
disp(d);
exts1={};
exts0={};
exts0205={};
exts0231={};
exts0951={};
for i=0:(Nz-1)
    disp(i);
    wpname0 = strcat('wpair',num2str(i));
    model.geom('geom1').create(wpname0, 'WorkPlane');
    wp0 = model.geom('geom1').feature(wpname0);
    wp0.set('unite', true);
    wp0.set('quickz', strcat(num2str(dist*i),'[nm]'));
    
    model.geom('geom1').run(wpname0);

    wpname1 = strcat('wpmat',num2str(i));
    model.geom('geom1').create(wpname1, 'WorkPlane');
    wp1 = model.geom('geom1').feature(wpname1);
    wp1.set('unite', true);
    wp1.set('quickz', strcat(num2str(dist*i),'[nm]'));

    model.geom('geom1').run(wpname1);

    wpname0205 = strcat('wp0205',num2str(i));
    model.geom('geom1').create(wpname0205, 'WorkPlane');
    wp0205 = model.geom('geom1').feature(wpname0205);
    wp0205.set('unite', true);
    wp0205.set('quickz', strcat(num2str(dist*i),'[nm]'));

    model.geom('geom1').run(wpname0205);

    wpname0231 = strcat('wp0231',num2str(i));
    model.geom('geom1').create(wpname0231, 'WorkPlane');
    wp0231 = model.geom('geom1').feature(wpname0231);
    wp0231.set('unite', true);
    wp0231.set('quickz', strcat(num2str(dist*i),'[nm]'));

    model.geom('geom1').run(wpname0231);

    wpname0951 = strcat('wp0951',num2str(i));
    model.geom('geom1').create(wpname0951, 'WorkPlane');
    wp0951 = model.geom('geom1').feature(wpname0951);
    wp0951.set('unite', true);
    wp0951.set('quickz', strcat(num2str(dist*i),'[nm]'));

    model.geom('geom1').run(wpname0951);

    rects1 = {};
    rects0 = {};
    rects0205={};
    rects0231={};
    rects0951={};
    for j=1:Nx*Ny
        %disp(bitmask(i*Nx*Ny+j));
        if bitmask(3*(i*Nx*Ny+j))==1
            rectname = strcat('rmat', num2str(Nx*Ny*i+j));
            rects1{length(rects1)+1} = rectname;
            wp1.geom().create(rectname, 'Rectangle');
            wp1.geom().feature(rectname).set('size', {d,d});
            x = strcat(num2str(mod(j-1,Nx)*dist),'[nm]');
            y = strcat(num2str(ceil(j/Nx)*dist),'[nm]');
           % disp(x);
            wp1.geom().feature(rectname).set('pos', {x,y});
            wp1.geom().feature(rectname).set('base', 'corner');

        elseif bitmask(3*(i*Nx*Ny+j))==0
            % disp('entered 0!');
            rectname = strcat('rair', num2str(Nx*Ny*i+j));
            rects0{length(rects0)+1} = rectname;
            wp0.geom().create(rectname, 'Rectangle');
            wp0.geom().feature(rectname).set('size', {d,d});
            x = strcat(num2str(mod(j-1,Nx)*dist),'[nm]');
            y = strcat(num2str(ceil(j/Nx)*dist),'[nm]');
           % disp(x);
            wp0.geom().feature(rectname).set('pos', {x,y});
            wp0.geom().feature(rectname).set('base', 'corner');

        elseif bitmask(3*(i*Nx*Ny+j))==0.205
            disp('entered 0.205!');
            rectname = strcat('r0205', num2str(Nx*Ny*i+j));
            rects0205{length(rects0205)+1} = rectname;
            wp0205.geom().create(rectname, 'Rectangle');
            wp0205.geom().feature(rectname).set('size', {d,d});
            x = strcat(num2str(mod(j-1,Nx)*dist),'[nm]');
            y = strcat(num2str(ceil(j/Nx)*dist),'[nm]');
           % disp(x);
            wp0205.geom().feature(rectname).set('pos', {x,y});
            wp0205.geom().feature(rectname).set('base', 'corner');

        elseif bitmask(3*(i*Nx*Ny+j))==0.231
            disp('entered 0.231!');
            rectname = strcat('r0231', num2str(Nx*Ny*i+j));
            rects0231{length(rects0231)+1} = rectname;
            wp0231.geom().create(rectname, 'Rectangle');
            wp0231.geom().feature(rectname).set('size', {d,d});
            x = strcat(num2str(mod(j-1,Nx)*dist),'[nm]');
            y = strcat(num2str(ceil(j/Nx)*dist),'[nm]');
           % disp(x);
            wp0231.geom().feature(rectname).set('pos', {x,y});
            wp0231.geom().feature(rectname).set('base', 'corner');

        elseif bitmask(3*(i*Nx*Ny+j))==0.951
            disp('entered 0.951!');
            rectname = strcat('rair', num2str(Nx*Ny*i+j));
            rects0951{length(rects0951)+1} = rectname;
            wp0951.geom().create(rectname, 'Rectangle');
            wp0951.geom().feature(rectname).set('size', {d,d});
            x = strcat(num2str(mod(j-1,Nx)*dist),'[nm]');
            y = strcat(num2str(ceil(j/Nx)*dist),'[nm]');
           % disp(x);
            wp0951.geom().feature(rectname).set('pos', {x,y});
            wp0951.geom().feature(rectname).set('base', 'corner');
        else

        end
    end
    uniname1 = strcat('unil1',num2str(i));
    wp1.geom().create(uniname1, 'Union');
    wp1.geom().feature(uniname1).selection('input').set(rects1);
    wp1.geom().feature(uniname1).set('intbnd', false);
    extname1 = strcat('extmat',num2str(i));
    exts1{length(exts1)+1}=extname1;
    model.geom('geom1').create(extname1, 'Extrude');
    model.geom('geom1').feature(extname1).set('workplane', wpname1);
    disp(i);
    model.geom('geom1').feature(extname1).selection('input').set(wpname1);
    model.geom('geom1').feature(extname1).setIndex('distance', d, 0);

    uniname0 = strcat('unil0',num2str(i));
    wp0.geom().create(uniname0, 'Union');
    wp0.geom().feature(uniname0).selection('input').set(rects0);
    wp0.geom().feature(uniname0).set('intbnd', false);
    extname0 = strcat('extair',num2str(i));
    exts0{length(exts0)+1}=extname0;
    model.geom('geom1').create(extname0, 'Extrude');
    model.geom('geom1').feature(extname0).set('workplane', wpname0);
    model.geom('geom1').feature(extname0).selection('input').set(wpname0);
    model.geom('geom1').feature(extname0).setIndex('distance', d, 0);

    uniname0231 = strcat('unil0231',num2str(i));
    disp('starting with 0.231');
    wp0231.geom().create(uniname0231, 'Union');
    wp0231.geom().feature(uniname0231).selection('input').set(rects0231);
    wp0231.geom().feature(uniname0231).set('intbnd', false);
    extname0231 = strcat('ext0231',num2str(i));
    exts0231{length(exts0231)+1}=extname0231;
    model.geom('geom1').create(extname0231, 'Extrude');
    model.geom('geom1').feature(extname0231).set('workplane', wpname0231);
    model.geom('geom1').feature(extname0231).selection('input').set(wpname0231);
    model.geom('geom1').feature(extname0231).setIndex('distance', d, 0);

    uniname0205 = strcat('unil0205',num2str(i));
    disp('starting 0.205')
    wp0205.geom().create(uniname0205, 'Union');
    wp0205.geom().feature(uniname0205).selection('input').set(rects0205);
    wp0205.geom().feature(uniname0205).set('intbnd', false);
    extname0205 = strcat('ext0205',num2str(i));
    exts0205{length(exts0205)+1}=extname0205;
    model.geom('geom1').create(extname0205, 'Extrude');
    model.geom('geom1').feature(extname0205).set('workplane', wpname0205);
    model.geom('geom1').feature(extname0205).selection('input').set(wpname0205);
    model.geom('geom1').feature(extname0205).setIndex('distance', d, 0);

    disp('finished with 0.205')

    uniname0951 = strcat('unil0951',num2str(i));
    wp0951.geom().create(uniname0951, 'Union');
    wp0951.geom().feature(uniname0951).selection('input').set(rects0951);
    wp0951.geom().feature(uniname0951).set('intbnd', false);
    extname0951 = strcat('ext0951',num2str(i));
    exts0951{length(exts0951)+1}=extname0951;
    model.geom('geom1').create(extname0951, 'Extrude');
    model.geom('geom1').feature(extname0951).set('workplane', wpname0951);
    model.geom('geom1').feature(extname0951).selection('input').set(wpname0951);
    model.geom('geom1').feature(extname0951).setIndex('distance', d, 0);
end
model.geom('geom1').create('uni1', 'Union');
model.geom('geom1').feature('uni1').set('intbnd', 'off');
model.geom('geom1').feature('uni1').selection('input').set(exts1);

model.geom('geom1').create('uni0', 'Union');
model.geom('geom1').feature('uni0').set('intbnd', 'off');
model.geom('geom1').feature('uni0').selection('input').set(exts0);

model.geom('geom1').create('uni0205', 'Union');
model.geom('geom1').feature('uni0205').set('intbnd', 'off');
model.geom('geom1').feature('uni0205').selection('input').set(exts0205);

model.geom('geom1').create('uni0231', 'Union');
model.geom('geom1').feature('uni0231').set('intbnd', 'off');
model.geom('geom1').feature('uni0231').selection('input').set(exts0231);

model.geom('geom1').create('uni0951', 'Union');
model.geom('geom1').feature('uni0951').set('intbnd', 'off');
model.geom('geom1').feature('uni0951').selection('input').set(exts0951);

model.geom.run();


model.save('E:\Topo-DDA-forWin64\Topo-DDA-forWin64\DDA_geo_cylinder_it300_piecewise_aboslute_unbinarized.mph')
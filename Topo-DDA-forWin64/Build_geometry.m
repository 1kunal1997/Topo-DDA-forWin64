name = "DDA_template.mph";
model = mphload(name);
data = importdata('cylinder_it300_lam542_sym_epsilon_0.1_penaltytype_piecewise_absolute_0.0to0.5_repeat\CoreStructure\CoreStructure299Rounded.txt');
bitmask = round(data);
Nx = 22;
Ny = 22;
Nz = 10;
dist = 15;
d=strcat(num2str(dist),'[nm]');
disp(d);
exts1={};
exts0={};
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
    rects1 = {};
    rects0 = {};
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
        else
            rectname = strcat('rair', num2str(Nx*Ny*i+j));
            rects0{length(rects0)+1} = rectname;
            wp0.geom().create(rectname, 'Rectangle');
            wp0.geom().feature(rectname).set('size', {d,d});
            x = strcat(num2str(mod(j-1,Nx)*dist),'[nm]');
            y = strcat(num2str(ceil(j/Nx)*dist),'[nm]');
           % disp(x);
            wp0.geom().feature(rectname).set('pos', {x,y});
            wp0.geom().feature(rectname).set('base', 'corner');
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
end
model.geom('geom1').create('uni1', 'Union');
model.geom('geom1').feature('uni1').set('intbnd', 'off');
model.geom('geom1').feature('uni1').selection('input').set(exts1);

model.geom('geom1').create('uni0', 'Union');
model.geom('geom1').feature('uni0').set('intbnd', 'off');
model.geom('geom1').feature('uni0').selection('input').set(exts0);

model.geom.run();


model.save('E:\Topo-DDA-forWin64\Topo-DDA-forWin64\DDA_geo_cylinder_it300_piecewise_aboslute_unbinarized.mph')
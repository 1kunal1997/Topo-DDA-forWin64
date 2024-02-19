name = 'DDA_template.mph';
model = mphload(name);
data = importdata('thick400-diel2d5-phi0theta0-lam500-size2000-focus450-2D_CoreStructure99.txt');
bitmask = data(4:end,1);
Nx = data(1,1);
Ny = data(2,1);
Nz = data(3,1);
dist = (data(1,3)-data(1,2))/(Nx-1);
d=strcat(num2str(dist),'[nm]');
exts={};
for i=1:9
    disp(i);
    wpname = strcat('wp',num2str(i));
    model.geom('geom1').create(wpname, 'WorkPlane');
    wp = model.geom('geom1').feature(wpname);
    wp.set('unite', true);
    wp.set('quickz', strcat(num2str(dist*i),'[nm]'));
    rects = {};
    for j=1:Nx*Ny
        if bitmask(i*Nx*Ny+j)==1
            rectname = strcat('r', num2str(Nx*Ny*i+j));
            rects{length(rects)+1} = rectname;
            wp.geom().create(rectname, 'Rectangle');
            wp.geom().feature(rectname).set('size', {d,d});
            x = strcat(num2str(mod(j,Nx)*dist),'[nm]');
            y = strcat(num2str(ceil(j/Nx)*dist),'[nm]');
            wp.geom().feature(rectname).set('pos', {x,y});
            wp.geom().feature(rectname).set('base', 'center');
        end
    end
    uniname = strcat('unil',num2str(i));
    wp.geom().create(uniname, 'Union');
    wp.geom().feature(uniname).selection('input').set(rects);
    wp.geom().feature(uniname).set('intbnd', false);
    extname = strcat('ext',num2str(i));
    exts{length(exts)+1}=extname;
    model.geom('geom1').create(extname, 'Extrude');
    model.geom('geom1').feature(extname).selection('input').set(wpname);
    model.geom('geom1').feature(extname).setIndex('distance', d, 0);
end
model.geom('geom1').create('uni1', 'Union');
model.geom('geom1').feature('uni1').set('intbnd', 'off');
model.geom('geom1').feature('uni1').selection('input').set(exts);
model.geom.run();


model.save('C:\Users\yz128\Documents\DDA\COMSOL verification\DDA_geo.mph')
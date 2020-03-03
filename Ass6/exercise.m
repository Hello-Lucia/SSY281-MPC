A = [0.5 0; 1 -0.5];


model = LTISystem('A',A);

model.x.min = [-10; -10];
model.x.max = [10; 10];

x = Polyhedron('lb', model.x.min, 'ub', model.x.max);

Xpre = model.reachableSet('X', x, 'direction', 'backward', 'N', 1);


figure(1)
plot([Xpre x])
y


%%
Xreach = model.reachableSet('X', x, 'direction', 'forward', 'N', 1);

figure(2)
plot([x Xreach])
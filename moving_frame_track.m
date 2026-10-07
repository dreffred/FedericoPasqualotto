%% MOVING_FRAME_TRACK
% Animates two orthonormal reference systems:
%   - a FIXED frame  {O; i, j, k}  at the origin
%   - a MOVING frame {P; T, N, B}  whose origin P runs in a loop along a 3D
%     parametric track r(t). The first unit vector T is always parallel
%     (tangent) to the track.
%
% ------------------------------------------------------------------------
% Theory (Frenet-Serret frame)
%   Velocity along the track:          r'(t)  = dr/dt
%   Unit tangent (parallel to track):  T = r' / |r'|
%   The plane spanned by r' and r'' is the osculating plane, so its normal
%   gives the unit binormal:           B = (r' x r'') / |r' x r''|
%   The principal normal completes a right-handed triad:
%                                      N = B x T
%   Curvature:                         kappa = |r' x r''| / |r'|^3   [1/m]
%
%   Since T, N, B are unit and mutually orthogonal, the matrix whose columns
%   are their components in the fixed frame
%                                      R = [T N B]
%   is a rotation matrix (R'*R = I, det R = +1): it maps moving-frame
%   components to fixed-frame components, v_fixed = R * v_moving.
%
%   Where kappa -> 0 (straight segments, inflection points) r' x r'' -> 0
%   and N, B are undefined. The script then keeps the previous N, projected
%   onto the plane normal to T, so the frame stays continuous.
%
% Alternative frame ('rmf'): rotation-minimizing (Bishop) frame, computed
% with the double-reflection method (Wang et al., ACM TOG 2008). T is still
% tangent, but N does not spin around T with the torsion. It is defined
% everywhere (also on straight lines) and is what CAD sweeps / camera rigs
% often use. On a closed track it generally does NOT return to its start
% orientation after one lap (holonomy), so it keeps propagating lap by lap.
%
% Controls: SPACE = pause/resume, Q (or close the window) = quit.
% Compatible with MATLAB (R2014b+) and GNU Octave.
% ------------------------------------------------------------------------

clear; clc; close all;

%% ---------------- USER SETTINGS ----------------
% Parametric track r(t) = [x(t); y(t); z(t)]   [m], t in [tStart, tEnd]
% Must return a 3x1 column for scalar t.
% Default: trefoil knot (closed curve, curvature never zero).
r = @(t) [ sin(t) + 2*sin(2*t);
           cos(t) - 2*cos(2*t);
          -sin(3*t) ];
tStart = 0;  tEnd = 2*pi;

% Other examples (uncomment one block):
% r = @(t) [2*cos(t); 2*sin(t); 0.3*t];  tStart = 0; tEnd = 4*pi;      % helix (open: restarts at the beginning)
% r = @(t) [1+cos(t); sin(t); 2*sin(t/2)]; tStart = 0; tEnd = 4*pi;     % Viviani curve
% r = @(t) [(2+cos(3*t))*cos(2*t); (2+cos(3*t))*sin(2*t); sin(3*t)]; tStart = 0; tEnd = 2*pi;  % (2,3) torus knot
% r = @(t) [3*cos(t); 2*sin(t); 0.8*sin(2*t)]; tStart = 0; tEnd = 2*pi; % saddle-shaped loop

frameType = 'frenet';  % 'frenet' or 'rmf'
nSteps    = 400;       % animation steps per lap [-]
Lfix      = 1.5;       % drawn length of the fixed unit vectors  [m] (display only)
Lmov      = 1.0;       % drawn length of the moving unit vectors [m] (display only)
pauseTime = 0.02;      % delay between frames [s]
kappaTol  = 1e-6;      % below this curvature [1/m] Frenet N,B are undefined

%% ---------------- PRECOMPUTATION ----------------
% Parameter grid: last point dropped so a closed curve is not repeated
tt = linspace(tStart, tEnd, nSteps + 1);
tt(end) = [];

% Track points (for drawing the path and setting the axes limits)
pts = zeros(3, nSteps);
for k = 1:nSteps
    pts(:, k) = r(tt(k));
end

% Numerical derivatives (central differences, error O(h^2)), so any r(t)
% works without differentiating it by hand.
hd  = 1e-4 * (tEnd - tStart);
dr  = @(t) (r(t + hd) - r(t - hd)) / (2*hd);
ddr = @(t) (r(t + hd) - 2*r(t) + r(t - hd)) / hd^2;

% Initial normal: Frenet N at tStart if defined, otherwise any unit vector
% perpendicular to T (fixed axis least aligned with T, projected).
rp = dr(tt(1));  rpp = ddr(tt(1));
T0 = rp / norm(rp);
c  = cross(rp, rpp);
if norm(c) / norm(rp)^3 > kappaTol
    Nprev = cross(c / norm(c), T0);
else
    E = eye(3);
    [~, idx] = min(abs(T0));
    Nprev = E(:, idx) - dot(E(:, idx), T0) * T0;
    Nprev = Nprev / norm(Nprev);
end
Tprev = T0;
xprev = pts(:, 1);

%% ---------------- FIGURE ----------------
fig = figure('Name', 'Fixed and moving frames', 'Color', 'w', ...
             'NumberTitle', 'off');
setappdata(fig, 'key', '');
set(fig, 'KeyPressFcn', @(src, ev) setappdata(src, 'key', ev.Key));
hold on; grid on; axis equal; view(40, 25);
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');

% Axis limits: track + fixed frame + room for the moving arrows
lo = min([pts, zeros(3,1)], [], 2) - max(Lmov, Lfix);
hi = max([pts, zeros(3,1)], [], 2) + max(Lmov, Lfix);
axis([lo(1) hi(1) lo(2) hi(2) lo(3) hi(3)]);

% Track
plot3(pts(1,:), pts(2,:), pts(3,:), '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

% Fixed frame {O; i, j, k}  (dashed, darker colours)
colFix = [0.6 0 0; 0 0.5 0; 0 0 0.6];
Efix   = eye(3);
fixNames = {'i', 'j', 'k'};
for a = 1:3
    quiver3(0, 0, 0, Lfix*Efix(1,a), Lfix*Efix(2,a), Lfix*Efix(3,a), 0, ...
            'Color', colFix(a,:), 'LineWidth', 2, 'LineStyle', '--', ...
            'MaxHeadSize', 0.4);
    text(1.1*Lfix*Efix(1,a), 1.1*Lfix*Efix(2,a), 1.1*Lfix*Efix(3,a), ...
         fixNames{a}, 'Color', colFix(a,:), 'FontWeight', 'bold');
end
plot3(0, 0, 0, 'ko', 'MarkerFaceColor', 'k');
text(0, 0, -0.2*Lfix, 'O');

% Moving frame {P; T, N, B} (solid, bright colours) - handles updated in loop
colMov = [1 0 0; 0 0.8 0; 0 0.4 1];
movNames = {'T', 'N', 'B'};
hMov = zeros(1, 3);
hTxt = zeros(1, 3);
for a = 1:3
    hMov(a) = quiver3(0, 0, 0, 0, 0, 0, 0, 'Color', colMov(a,:), ...
                      'LineWidth', 2.5, 'MaxHeadSize', 0.5);
    hTxt(a) = text(0, 0, 0, movNames{a}, 'Color', colMov(a,:), ...
                   'FontWeight', 'bold');
end
hP = plot3(0, 0, 0, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 7);
hTitle = title('');

%% ---------------- ANIMATION LOOP ----------------
k = 0;
paused = false;
while ishghandle(fig)
    % Keyboard controls
    key = getappdata(fig, 'key');
    if ~isempty(key)
        setappdata(fig, 'key', '');
        if strcmpi(key, 'space'), paused = ~paused; end
        if strcmpi(key, 'q'),     close(fig); break;  end
    end
    if paused
        pause(0.05);
        continue;
    end

    k = mod(k, nSteps) + 1;          % loop index: 1..nSteps, 1..nSteps, ...
    t = tt(k);
    x = pts(:, k);                   % origin P of the moving frame [m]

    rp  = dr(t);                     % r'  [m per unit t]
    rpp = ddr(t);                    % r'' [m per unit t^2]
    T   = rp / norm(rp);             % unit tangent -> parallel to the track
    c   = cross(rp, rpp);
    kappa = norm(c) / norm(rp)^3;    % curvature [1/m]

    switch lower(frameType)
        case 'frenet'
            if kappa > kappaTol
                B = c / norm(c);
                N = cross(B, T);
                % Note: N flips sign when the curve passes an inflection;
                % that is a property of the Frenet frame, not a bug.
            else
                % Undefined Frenet normal: carry the previous one over
                N = Nprev - dot(Nprev, T) * T;
                N = N / norm(N);
                B = cross(T, N);
            end

        case 'rmf'
            % Double reflection: reflect (Nprev, Tprev) across the plane
            % bisecting xprev->x, then across the plane bisecting TL->T.
            v1 = x - xprev;   c1 = dot(v1, v1);
            if c1 > eps
                NL = Nprev - (2/c1) * dot(v1, Nprev) * v1;
                TL = Tprev - (2/c1) * dot(v1, Tprev) * v1;
            else
                NL = Nprev;  TL = Tprev;
            end
            v2 = T - TL;      c2 = dot(v2, v2);
            if c2 > eps
                N = NL - (2/c2) * dot(v2, NL) * v2;
            else
                N = NL;
            end
            % Re-orthonormalise to stop round-off drift over many laps
            N = N - dot(N, T) * T;
            N = N / norm(N);
            B = cross(T, N);

        otherwise
            error('frameType must be ''frenet'' or ''rmf''.');
    end
    Nprev = N;  Tprev = T;  xprev = x;

    % Rotation matrix moving -> fixed and orthonormality check
    R = [T, N, B];
    orthoErr = norm(R' * R - eye(3));

    % Update graphics (window may have been closed during pause)
    if ~ishghandle(fig), break; end
    for a = 1:3
        set(hMov(a), 'XData', x(1), 'YData', x(2), 'ZData', x(3), ...
                     'UData', Lmov*R(1,a), 'VData', Lmov*R(2,a), ...
                     'WData', Lmov*R(3,a));
        set(hTxt(a), 'Position', (x + 1.15*Lmov*R(:,a))');
    end
    set(hP, 'XData', x(1), 'YData', x(2), 'ZData', x(3));
    set(hTitle, 'String', { ...
        sprintf('%s frame   t = %.3f   \\kappa = %.3f 1/m   |R^TR - I| = %.1e   det R = %+.3f', ...
                upper(frameType), t, kappa, orthoErr, det(R)), ...
        sprintf('T = [%+.2f %+.2f %+.2f]   N = [%+.2f %+.2f %+.2f]   B = [%+.2f %+.2f %+.2f]', ...
                T, N, B)});

    drawnow;
    pause(pauseTime);
end

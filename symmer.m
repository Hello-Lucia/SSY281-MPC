function print = symmer(stuff)
    clf;
    print = latex(stuff);
    
    % p
    print = regexprep(print, 'ddp(?!hi)', '\\ddot\{p\}');
    print = regexprep(print, 'dp(?!hi)', '\\dot\{p\}');
    
    %theta
    print = strrep(print, 'ddtheta', '\ddot{\theta}');
    print = strrep(print, 'dtheta', '\dot{\theta}');
    
    %alpha
    print = strrep(print, 'ddalpha', '\ddot{\alpha}');
    print = strrep(print, 'dalpha', '\dot{\alpha}');
    
    %phi
    print = strrep(print, 'ddphi', '\ddot{\phi}');
    print = strrep(print, 'dphi', '\dot{\phi}');
    
    % Small phi
    print = strrep(print, '\phi', '\varphi');
    
    
    math_print = strcat("$", print, "$");
    
    axis off;
    text(0.5, 0.5, math_print, ...
        'HorizontalAlignment', 'Center', ...
        'VerticalAlignment', 'Middle', ...
        'Interpreter', 'latex', ...
        'FontSize', 20);
end
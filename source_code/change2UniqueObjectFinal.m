 function [struttura_nuova] = change2UniqueObjectFinal(struttura)

    N_objects = length(struttura.aggregate);
    struttura_nuova.fv.faces = [];
    struttura_nuova.fv.vertices = [];
    
    for i=1:N_objects
        
        dim = size(struttura.aggregate(i).s_Object.fv.faces);
        n_faces(i) = dim(1);
        struttura_nuova.fv.vertices = [struttura_nuova.fv.vertices; ...
                                       struttura.aggregate(i).s_Object.fv.vertices];
                                  
    end
    
    vettore_faces = [1:1:sum(n_faces)*3];
    struttura_nuova.fv.faces = reshape(vettore_faces, [3,sum(n_faces)])';

end
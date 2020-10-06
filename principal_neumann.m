% =====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions de
% Neumann sur le maillage nom_maillage.msh
%
% | -\Delta u + u= f,   dans \Omega
% |         du/dn = 0,   sur le bord
%
% =====================================================


% lecture du maillage et affichage
% ---------------------------------
nom_maillage = 'geomCarre.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% affichemaillage(nom_maillage);

matM_elem([0, 0], [0, 1], [1, 1]);

% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
LL = zeros(Nbpt,1);     % vecteur second membre

% boucle sur les triangles
% ------------------------
for l = 1:Nbtri
  % Coordonnees des sommets du triangles
  % A COMPLETER
  S1 = Coorneu(Numtri(l, 1),:);
  S2 = Coorneu(Numtri(l, 2),:);
  S3 = Coorneu(Numtri(l, 3),:);
  % calcul des matrices elementaires du triangle l 
  
  Kel = matK_elem(S1, S2, S3);
           
  Mel = matM_elem(S1, S2, S3);
    
  % On fait l'assemmblage de la matrice de masse et de rigidité
  for i = 1:3
    I = Numtri(l, i);
    for j = 1:3
      J = Numtri(l, j);
      MM(I,J) += Mel(i,j);
      KK(I,J) += Kel(i,j);
    endfor
  endfor
endfor

% Calcul du second membre L
% -------------------------
	% A COMPLETER
	% utiliser la routine f.m
  
FF = zeros(Nbpt,1);
for I = 1:Nbpt
  FF(I) = f(Coorneu(I,1), Coorneu(I,2));
endfor

LL = MM * FF;

% inversion
% ----------
UU = (MM+KK)\LL;

% visualisation
% -------------
affiche(UU, Numtri, Coorneu, sprintf('Neumann - %s', nom_maillage));

validation = 'non';
% validation
% ----------
if strcmp(validation,'oui')
UU_exact = cos(pi*Coorneu(:,1)).*cos(2*pi*Coorneu(:,2));
% Calcul de l erreur L2
% A COMPLETER
% Calcul de l erreur H1
% A COMPLETER
% attention de bien changer le terme source (dans FF)
endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020
% André Bastos
% Code to re-index PSD/CSD matricies by depth across sessions
%24 contqacts at spacing of 100um
coordinate_axis = -3200:100:3200;
dist_to_top= sink-1;
dist_to_bottom = 24-sink+1;
top_coord = dist_to_top*-100;
bottom_coord = dist_to_bottom*100;
new_coords = top_coord:100:bottom_coord; %cop into ElectrodeInfo
indx1 = find(coordinate_axis== top_coord)
indx2 = find(coordinate_axis== bottom_coord)
master_matrix_POW(n,indx1:indx2,:) = this_session_POW;
master_matrix_CSD(n,indx1:indx2,:) = this_session_CSD;
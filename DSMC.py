import numpy as np
import random


k_B  = 1.0
mu_B = 1.0
gs   = 1.0
mass = 1.0

N_ATOMS = 10
TEMP = 20. * 10 ** -6

max_grid_width = 1.0
CELLS_PER_AXIS = 2
NUM_OF_CELLS = CELLS_PER_AXIS ** 3
N_TH = 5

def get_gaussian_point( mean,
						std ):

	rndpt    = np.empty( 3 )
	rndpt[0] = random.gauss( mean, 
			  			     std )
	rndpt[1] = random.gauss( mean, 
			  			     std )
	rndpt[2] = random.gauss( mean, 
			  			     std )

	return rndpt

def get_random_velocity( temp ):

	V = np.sqrt( k_B * temp / mass )

	vel = get_gaussian_point( 0.,
		    				  V )

	return vel

def select_pos_in_thermal_dist( temp ):

	R = np.sqrt( k_B * temp / mass )

	pos = get_gaussian_point( 0.,
		    				  R )

	# no_atom_selected = True

	# while no_atom_selected:
	# 	r = get_gaussian_point( 0.,
	# 		                    1. )
	# 	r = r * max_grid_width / 3.

	# 	magB = norm( B( r ) )
	# 	U = 0.5 * ( magB - B0 ) * gs * mu_B
	# 	Pr = np.exp( -U / k_B / temp )

	# 	if random.random() < Pr:
	# 		pos = r
	# 		no_atom_selected = False

	return pos

def generate_initial_dist( number_of_atoms,
			   			   initial_temp ):

	pos     = [ ]
	vel     = [ ]
	atom_id = [ ]

	for atom in range(number_of_atoms):
		pos.append( select_pos_in_thermal_dist( initial_temp ) )
		vel.append( get_random_velocity( initial_temp ) )
		atom_id.append( number_of_atoms - atom - 1 )

	return np.array(pos), np.array(vel), np.array(atom_id)

def sort( pos,
	      atom_id ):

	

	return

def collide_atoms( pos,
	               vel,
	               cell_min,
	               cell_width,
	               dt,
	               sig_vr_max,
	               cell_start_end,
	               cell_id,
	               atom_id,
	               number_of_collisions,
	               num_atoms,
	               alpha,
	               atom_count,
	               num_cells,
	               level ):

	for cell in range(num_cells):
		print( cell ) 
		new_level = level + 1
		
		l_cell_start_end = cell_start_end[cell]
		l_cell_id = cell_id[l_cell_start_end[0]:l_cell_start_end[1]+1]
		l_atom_id = atom_id[l_cell_start_end[0]:l_cell_start_end[1]+1]
		num_atoms_in_cell  = get_number_of_atoms( l_cell_start_end )
		
		sub_cell_width = cell_width / CELLS_PER_AXIS
		sub_cell_min   = cell_min + sub_cell_width * get_sub_cell_index( cell,
																         CELLS_PER_AXIS)
		sub_cell_max   = sub_cell_min + cell_width

		if num_atoms_in_cell > N_TH:
			# print( "cell{0}-{1}, sub_cell_min = {2}".format(cell,new_level,sub_cell_min) )
			set_cell_id( pos,
				 		 l_cell_id,
				 		 l_atom_id,
				 		 sub_cell_min,
				 		 sub_cell_width,
				 		 CELLS_PER_AXIS,
				 		 num_atoms_in_cell )

			ind = np.lexsort((l_atom_id, l_cell_id))
			l_atom_id = l_atom_id[ind]
			l_cell_id = l_cell_id[ind]
			print( l_cell_id )

			# collide_atoms( pos,
	  #              	   vel,
	  #                  sub_cell_min,
	  #                  sub_cell_width,
	  #                  dt,
	  #                  sig_vr_max,
	  #                  cell_start_end,
	  #                  cell_id,
	  #                  atom_id,
	  #                  number_of_collisions,
	  #                  num_atoms,
	  #                  alpha,
	  #                  atom_count,
	  #                  NUM_OF_CELLS,
	  #                  new_level )
		else:
			return
			# do_a_collision()

	return

def set_cell_id( pos,
				 cell_id,
				 atom_id,
				 cell_min,
				 cell_width,
				 cells_per_axis,
				 number_of_atoms ):

	for atom in range( number_of_atoms ):

		cell_index = get_cell_index( pos[atom_id[atom]],
								     cell_min,
								     cell_width )

		cell_id[atom] = get_cell_id( cell_index,
								     cells_per_axis )

def get_cell_index( pos,
				    cell_min,
				    cell_width ):
	index = np.zeros(3,).astype(np.int)

	index[0] = ( pos[0] - cell_min[0] ) / cell_width[0]
	index[1] = ( pos[1] - cell_min[1] ) / cell_width[1]
	index[2] = ( pos[2] - cell_min[2] ) / cell_width[2]
	
	return index

def get_cell_id( index,
				 cells_per_axis ):
	cell_id = 0
	
	if( index[0] > -1 and index[0] < cells_per_axis and
		index[1] > -1 and index[1] < cells_per_axis and
		index[2] > -1 and index[2] < cells_per_axis ):
		cell_id = index[2]*cells_per_axis*cells_per_axis + index[1]*cells_per_axis + index[0]
	else:
		cell_id = cells_per_axis * cells_per_axis * cells_per_axis
	
	return cell_id

def get_sub_cell_index( cell_id,
					    cells_per_axis ):
	index = np.zeros(3,).astype(np.int)

	index[2] = cell_id / (cells_per_axis*cells_per_axis)
	index[1] = (cell_id - index[2]*(cells_per_axis*cells_per_axis )) / cells_per_axis
	index[0] = cell_id - index[2]*cells_per_axis*cells_per_axis - index[1]*cells_per_axis

	return index

def get_number_of_atoms( cell_start_end ):

	number_of_atoms  = cell_start_end[1] - cell_start_end[0] + 1
	
	if (cell_start_end[0] < 0 or cell_start_end[1] < 0):
		number_of_atoms = 0

	return number_of_atoms

def main():

	pos, vel, atom_id = generate_initial_dist( N_ATOMS,
			           		   				   TEMP )
	
	print( len(atom_id) )

	cell_id = np.zeros(N_ATOMS,).astype(np.int)
	set_cell_id( pos,
				 cell_id,
				 atom_id,
				 np.array([-5.,-5.,-5.]),
				 np.array([10.,10.,10.]),
				 CELLS_PER_AXIS,
				 N_ATOMS )

	ind = np.lexsort((atom_id, cell_id))
	atom_id = atom_id[ind]
	cell_id = cell_id[ind]

	collide_atoms( pos,
	               vel,
	               np.array([-5.,-5.,-5.]),
	               np.array([10.,10.,10.]),
	               None,
	               None,
	               np.array([[0,N_ATOMS-1],]),
	               cell_id,
	               atom_id,
	               0.,
	               N_ATOMS,
	               1.,
	               0,
	               1,
	               0 )

	# print( atom_id )

	return

if __name__ == "__main__":
    main()
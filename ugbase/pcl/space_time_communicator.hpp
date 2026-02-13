/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
 * Author: Arne Nägel
 *
 * This file is part of UG4.
 *
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__PCL__PCL_SPACE_TIME_COMMUNICATOR_H
#define __H__PCL__PCL_SPACE_TIME_COMMUNICATOR_H


#include "pcl/pcl_comm_world.h"

namespace pcl{

    class SpaceTimeCommunicator {
    public:

        //--------------------------------------------------------------------------------------------------------------

        SpaceTimeCommunicator() = default;
        virtual ~SpaceTimeCommunicator() = default;

        //--------------------------------------------------------------------------------------------------------------

        void split(int numTemporalProcesses) {

            int world_size, myid;
            MPI_Comm_size(PCL_COMM_WORLD, &world_size);
			if(world_size % numTemporalProcesses != 0 )
				UG_THROW("SpaceTimeCommunicator: Not enough processses for temporal spliting");
            

            GLOBAL = PCL_COMM_WORLD;

            globalsize_ = world_size;
            spatialsize_ = world_size / numTemporalProcesses;
            temporalsize_ = numTemporalProcesses;
            

            MPI_Comm_rank(GLOBAL, &myid);
            const int xcolor = myid / spatialsize_;
            const int tcolor = myid % spatialsize_;

            MPI_Comm_split(GLOBAL, xcolor, myid, &SPATIAL);
            MPI_Comm_split(GLOBAL, tcolor, myid, &TEMPORAL);
			
            if (verbose_) {
                std::cout << "World size after splitting is:\t" << world_size << std::endl;
                std::cout << "... with temporal world size:\t" << temporalsize_ << std::endl;
                std::cout << "... and spatial world size:\t" << spatialsize_ << std::endl << std::endl;
            }
            

            PCL_COMM_WORLD = SPATIAL; // replaces ugs world communicator with the communicator for spatial
        }

        void unsplit() {
            PCL_COMM_WORLD = GLOBAL; // reset the world communicator
			
			MPI_Comm_free(&SPATIAL); // free the spatial communicator
            SPATIAL = PCL_COMM_WORLD;
			
			MPI_Comm_free(&TEMPORAL);// free the temporal communicator
            TEMPORAL = PCL_COMM_WORLD;
        }

        int get_global_size() const {
            return globalsize_;
        }

        int get_temporal_size() const {
            return temporalsize_;
        }

        int get_spatial_size() const {
            return spatialsize_;
        }

        int get_temporal_rank() const {
            int rank = 0;
            MPI_Comm_rank(TEMPORAL, &rank);
            return rank;
        }

        int get_spatial_rank() const {
            int rank = 0;
            MPI_Comm_rank(SPATIAL, &rank);
            return rank;
        }

        int get_global_rank() const {
            int rank = 0;
            MPI_Comm_rank(GLOBAL, &rank);
            return rank;
        }

        void sleep(int microseconds) {
            usleep(microseconds);
        }

        //--------------------------------------------------------------------------------------------------------------
        MPI_Comm GLOBAL = PCL_COMM_WORLD;
        MPI_Comm TEMPORAL = PCL_COMM_WORLD;
        MPI_Comm SPATIAL = PCL_COMM_WORLD;

        bool verbose_ = true;
        int globalsize_ = 1;
        int temporalsize_ = 1;
        int spatialsize_ = 1;

        //--------------------------------------------------------------------------------------------------------------
    };
}
#endif

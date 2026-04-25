// Copyright 2019 Alexander Liniger

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//     http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#include "cost_test.h"
namespace mpcc{
int testCost(const PathToJson &path){
    Cost cost = Cost(path);
    //cost.setCosts(cost_param);
    ArcLengthSpline track = ArcLengthSpline(path);
    genRoundTrack(track);


    StateVector xk0_vec, xk1_vec, xk2_vec, xk3_vec, xk4_vec;
    xk0_vec <<  1, 0, 0, 2, 0, 0, 0      , 0, 0, 0;
    xk1_vec <<  0, 1, 0, 2, 0, 0, M_PI/2., 0, 0, 0;
    xk2_vec << -1, 0, 0, 2, 0, 0, M_PI   , 0, 0, 0;
    xk3_vec << -1.2, 0, 0, 2, 0, 0, M_PI , 0, 0, 0;
    xk4_vec << -1, 0, 0, 2, 0.5, 0, M_PI , 0, 0, 0;
    InputVector uk0_vec, uk1_vec, uk2_vec, uk3_vec, uk4_vec;
    uk0_vec << 0, 0, 0;
    uk1_vec << 0, 0, 0;
    uk2_vec << 0, 0, 0;
    uk3_vec << 1, 1, 1;
    uk4_vec << 0, 0, 0;


    //zero contouring error tests
    // contouring cost is approximatly the zero order term of the cost given the error is zero
    // 0.5 xQx + 0.5 uRu + qx + ru ~ 0.5 xQx
    // s = 0    -> right of circle
    CostMatrix cost_mat0 = cost.getCost(track,vectorToState(xk0_vec),1);

    if(std::fabs((0.5*xk0_vec.transpose()*cost_mat0.Q*xk0_vec + 0.5*uk0_vec.transpose()*cost_mat0.R*uk0_vec +  cost_mat0.q.transpose()*xk0_vec +  cost_mat0.r.transpose()*uk0_vec  + 0.5*xk0_vec.transpose()*cost_mat0.Q*xk0_vec).value()) >= 0.1){
        return 2;
    }
    // s = R*pi/2 -> top of circle
    CostMatrix cost_mat1 = cost.getCost(track,vectorToState(xk1_vec),1);


    if(std::fabs((0.5*xk1_vec.transpose()*cost_mat1.Q*xk1_vec + 0.5*uk1_vec.transpose()*cost_mat1.R*uk1_vec +  cost_mat1.q.transpose()*xk1_vec +  cost_mat1.r.transpose()*uk1_vec  + 0.5*xk1_vec.transpose()*cost_mat1.Q*xk1_vec).value()) >= 0.1) {
        return 3;
    }
    // s = R*pi -> left of circle
    CostMatrix cost_mat2 = cost.getCost(track,vectorToState(xk2_vec),1);

//    std::cout << std::fabs((0.5*xk2_vec.transpose()*Q*xk2_vec + 0.5*uk2_vec.transpose()*R*uk2_vec +  q.transpose()*xk2_vec +  r.transpose()*uk2_vec  + 0.5*xk2_vec.transpose()*Q*xk2_vec).value()) << std::endl;
    if(std::fabs((0.5*xk2_vec.transpose()*cost_mat2.Q*xk2_vec + 0.5*uk2_vec.transpose()*cost_mat2.R*uk2_vec +  cost_mat2.q.transpose()*xk2_vec +  cost_mat2.r.transpose()*uk2_vec  + 0.5*xk2_vec.transpose()*cost_mat2.Q*xk2_vec).value()) >= 0.1){
        return 4;
    }


    // Test offset from trajectory
    // given a point with error 0 does moving in X-Y change the error?
    if( (0.5*xk2_vec.transpose()*cost_mat2.Q*xk2_vec + 0.5*uk2_vec.transpose()*cost_mat2.R*uk2_vec +  cost_mat2.q.transpose()*xk2_vec +  cost_mat2.r.transpose()*uk2_vec).value() >=
        (0.5*xk3_vec.transpose()*cost_mat2.Q*xk3_vec + 0.5*uk2_vec.transpose()*cost_mat2.R*uk2_vec +  cost_mat2.q.transpose()*xk3_vec +  cost_mat2.r.transpose()*uk2_vec).value()){
        return 5;
    }
    xk3_vec(0) = 0.8;
    if( (0.5*xk2_vec.transpose()*cost_mat2.Q*xk2_vec + 0.5*uk2_vec.transpose()*cost_mat2.R*uk2_vec +  cost_mat2.q.transpose()*xk2_vec +  cost_mat2.r.transpose()*uk2_vec).value() >=
        (0.5*xk3_vec.transpose()*cost_mat2.Q*xk3_vec + 0.5*uk2_vec.transpose()*cost_mat2.R*uk2_vec +  cost_mat2.q.transpose()*xk3_vec +  cost_mat2.r.transpose()*uk2_vec).value()){
        return 6;
    }
    xk3_vec(0) = 1;
    xk3_vec(1) = 0.1;
    if( (0.5*xk2_vec.transpose()*cost_mat2.Q*xk2_vec + 0.5*uk2_vec.transpose()*cost_mat2.R*uk2_vec +  cost_mat2.q.transpose()*xk2_vec +  cost_mat2.r.transpose()*uk2_vec).value() >=
        (0.5*xk3_vec.transpose()*cost_mat2.Q*xk3_vec + 0.5*uk2_vec.transpose()*cost_mat2.R*uk2_vec +  cost_mat2.q.transpose()*xk3_vec +  cost_mat2.r.transpose()*uk2_vec).value()){
        return 7;
    }
    xk3_vec(1) = -0.1;
    if( (0.5*xk2_vec.transpose()*cost_mat2.Q*xk2_vec + 0.5*uk2_vec.transpose()*cost_mat2.R*uk2_vec +  cost_mat2.q.transpose()*xk2_vec +  cost_mat2.r.transpose()*uk2_vec).value() >=
        (0.5*xk3_vec.transpose()*cost_mat2.Q*xk3_vec + 0.5*uk2_vec.transpose()*cost_mat2.R*uk2_vec +  cost_mat2.q.transpose()*xk3_vec +  cost_mat2.r.transpose()*uk2_vec).value()){
        return 8;
    }

    // Test input increase
    // do increased inputs increase the cost?
    if( (0.5*xk2_vec.transpose()*cost_mat2.Q*xk2_vec + 0.5*uk2_vec.transpose()*cost_mat2.R*uk2_vec +  cost_mat2.q.transpose()*xk2_vec +  cost_mat2.r.transpose()*uk2_vec).value() >=
        (0.5*xk2_vec.transpose()*cost_mat2.Q*xk2_vec + 0.5*uk3_vec.transpose()*cost_mat2.R*uk3_vec +  cost_mat2.q.transpose()*xk2_vec +  cost_mat2.r.transpose()*uk3_vec).value()){
        return 9;
    }

    // Test Beta Cost
    // does an increased slip angle increase the cost?
    if( (0.5*xk2_vec.transpose()*cost_mat2.Q*xk2_vec + 0.5*uk2_vec.transpose()*cost_mat2.R*uk2_vec +  cost_mat2.q.transpose()*xk2_vec +  cost_mat2.r.transpose()*uk2_vec).value() >=
        (0.5*xk4_vec.transpose()*cost_mat2.Q*xk4_vec + 0.5*uk4_vec.transpose()*cost_mat2.R*uk4_vec +  cost_mat2.q.transpose()*xk4_vec +  cost_mat2.r.transpose()*uk4_vec).value()){
        return 10;
    }
    xk4_vec(4) = -0.5;
    if( (0.5*xk2_vec.transpose()*cost_mat2.Q*xk2_vec + 0.5*uk2_vec.transpose()*cost_mat2.R*uk2_vec +  cost_mat2.q.transpose()*xk2_vec +  cost_mat2.r.transpose()*uk2_vec).value() >=
        (0.5*xk4_vec.transpose()*cost_mat2.Q*xk4_vec + 0.5*uk4_vec.transpose()*cost_mat2.R*uk4_vec +  cost_mat2.q.transpose()*xk4_vec +  cost_mat2.r.transpose()*uk4_vec).value()){
        return 11;
    }

    return 1;
}
}
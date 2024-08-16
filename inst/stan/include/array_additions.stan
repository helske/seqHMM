  //
  // Helper functions for elementwise addition of integer arrays
  // By Rok Cesnovar at 
  // https://github.com/stan-dev/math/issues/2782#issuecomment-1178160026
  //
  array[] int add(array[] int x, array[] int y) {
    int x_size = size(x);
    array[x_size] int z;
    for (i in 1:x_size){
      z[i] = x[i] + y[i];
    }
    return z;
  }
  array[] int add(int x, array[] int y) {
    int y_size = size(y);
    array[y_size] int z;
    for (i in 1:y_size){
      z[i] = x + y[i];
    }
    return z;
  }
  array[] int add(array[] int x, int y) {
    int x_size = size(x);
    array[x_size] int z;
    for (i in 1:x_size){
      z[i] = x[i] + y;
    }
    return z;
  }

/*

Statistical methods

*/
module Statistics {

  import Sort.sort;

  /*
   Compute median of an array

   :arg arr: numeric array

   :returns: median of array
   */
  proc median(const arr: [?D] ?T): real
              where D.rank == 1 && isRealType(T) || isIntegralType(T) {
    var data = arr; // copy
    sort(data);
    var midIdx: int = (data.size-1) / 2;
    if data.size % 2 == 0 then return (data[midIdx] + data[midIdx + 1]) / 2.0;
    return data[midIdx];
  }

}
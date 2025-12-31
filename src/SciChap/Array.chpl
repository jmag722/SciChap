/*
  Array: Basic array operations.
 */
module Array {

  /*
  Computes difference between array elements

  :arg arr: input array

  :returns: Array of size `arr.size - 1`
  :rtype: [] T
  */
  proc diff(const ref arr: [?D] ?T): [0..<D.size-1] T
            where (D.rank == 1 && isNumeric(T)) {
    return arr[D.first+1..] - arr[..<D.last];
  }

  /*
  Determines if array is strictly increasing or decreasing

  :arg x: input array

  :returns: true if all elements strictly increasing or decreasing, else false
  :rtype: bool
  */
  proc isMonotonic(const ref x: [?D] ?T): bool
                   where (isNumeric(T) && !isComplex(T)) {
    if D.size < 2 then return true;
    const h: [0..<D.size - 1] real = diff(x);
    return (&& reduce (h > 0)) || (&& reduce (h < 0));
  }

  /*
  Creates linearly spaced array of elements

  :arg start: starting value
  :arg end: ending value
  :arg num: number of points in array
  :arg endpoint: if true, includes `end` in returned array, else excluded

  :returns: linearly spaced array of elements
  :rtype: [] T
  */
  proc linspace(in start: ?T, in end: T, in num: int,
                in endpoint: bool=true): [0..<num] T
                where (isNumeric(T) && !isIntegral(T)) {
    var arr: [0..<num] T;
    const divisor: int = if endpoint && num > 1 then num - 1 else num;
    const increment: T = (end - start) / divisor;
    forall i in arr.domain {
      arr[i] = i * increment + start;
    }
    return arr;
  }


  /*
   Create empty array

   :arg T: numeric array type
   :type T: type

   :returns: empty array
   :rtype: [] T
   */
  proc empty(type T): [1..0] T where isNumericType(T) {
    var a: [1..0] T;
    return a;
  }


  /*
   Converts a bounded range type into an array

   :arg rng: bounded input range
   :type rng: range

   :arg T: numeric type
   :type T: type

   :returns: range converted to array
   :rtype: [] T
   */
  proc arange(const rng: range(idxType=?, bounds=boundKind.both, strides=?),
              type T=int): [] T where isNumericType(T) {
    return [i in rng] i:T;
  }

  /*
   Finds all element indices in an array matching the specified value

   :arg x: input array
   :type x: [?D] ?T

   :arg value: specified value to find
   :type value: T

   :returns: indices where array values match input
   :rtype: [] int
   */
  proc findAll(const ref x: [?D] ?T, in value: T): [] int where D.rank == 1 {
    return [i in D] if x[i] == value then i;
  }

  /*
   Slice array with another array of booleans. Indexing an array with another
   array of booleans does not yield what you'd expect (it coerces boolean values
   to integers, then uses those as indices of the array).

   :arg x: input array
   :type x: [?D] ?T

   :arg mask: specified value to find
   :type mask: [D] bool

   :returns: array values where the mask was true
   :rtype: [] T
   */
  proc boolSlice(const ref x: [?D] ?T, const ref mask: [D] bool): [] T
                 where D.rank == 1 {
    return [i in D] if mask[i] then x[i];
  }

}
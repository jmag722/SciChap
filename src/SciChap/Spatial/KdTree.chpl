module KdTree {

  record kdTree {
    const dataDom: domain(2) = {1..0, 1..0};
    const data: [dataDom] real;

    proc init(const ref data: [?D] real): void where D.rank == 2 {
      this.dataDom = {0..#D.shape[0], 0..#D.shape[1]};
      this.data = data;
    }

  }

}
class Trace {
  constructor(
    mode,
    type,
    size,
    ids = [],
    x = [],
    y = [],
    cluster = [],
    z = []
  ) {
    this.text = ids;
    this.x = x;
    this.y = y;
    this.z = z;
    this.name = cluster;
    this.mode = mode;
    this.type = type;
    this.marker = { size: size };
  }
}

export { Trace };

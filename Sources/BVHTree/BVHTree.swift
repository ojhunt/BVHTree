//
//  BVHTree.swift
//  swiftrt
//
//  Created by Oliver Hunt on 4/2/21.
//

import Foundation
import VectorTypes

public protocol Intersectable {
  associatedtype PrimitiveType
  func intersect<RenderContextType: RenderContext>(context: RenderContextType,
                                                   ray: RenderContextType.RayType,
                                                   hitMode: HitMode,
                                                   min: RenderContextType.ValueType,
                                                   max: RenderContextType.ValueType) -> RenderContextType.CollisionType?
}
public protocol BVHElement {
  associatedtype PointType: VectorTypes.Point where PointType.VectorType : VectorTypes.Vector
  associatedtype LeafElementStorage: Intersectable
  associatedtype PrimitiveType
  static func buildLeafElements(_: [Self]) -> [LeafElementStorage]
  static func releaseLeafElements(_: [LeafElementStorage])
  var bounds: BoundingBox<PointType> { get }
}

public protocol IntConvertibleFloatingPoint: FloatingPoint, ExpressibleByFloatLiteral {
  var asInteger: Int { get }
  init(_: Float)
  init(_: Double)
}

extension Float: IntConvertibleFloatingPoint {
  public var asInteger: Int {
    return Int(self)
  }
}

extension Double: IntConvertibleFloatingPoint {
  public var asInteger: Int {
    return Int(self)
  }
}

public protocol Collision {
  associatedtype ElementType: BVHElement
  typealias ValueType = ElementType.PointType.ValueType
  var distance: ValueType { get }
  var target: ElementType.PrimitiveType { get }
  var uv: (ValueType, ValueType) { get }
  init(distance: ValueType,
       uv: (ValueType, ValueType),
       intersectionCount: Int,
       nodeCount: Int,
       target: ElementType.PrimitiveType);
}

extension BoundingBox where PointType: Point, PointType.ValueType: IntConvertibleFloatingPoint, PointType.VectorType: Vector {
  func offsetRatio(point: PointType, axis: PointType.AxisType) -> PointType.ValueType {
    let o = point - minBound;
    let scaleFactor = maxBound[axis] - minBound[axis];
    return o[axis] / scaleFactor
  }
  
  func surfaceArea() -> PointType.ValueType {
    let size = maxBound - minBound;
    var result = PointType.ValueType(1.0)
    for axis in PointType.AxisType.allAxes {
      result *= size[axis]
    }
    return result;
  }
  
  @inline(__always) func intersect<RayType: Ray>(_ ray: RayType, min near: RayType.ValueType, max far: RayType.ValueType) -> (min: RayType.ValueType, max: RayType.ValueType)?
  where RayType.PointType == PointType {
    var tmin = PointType.VectorType(repeating: near);
    var tmax = PointType.VectorType(repeating: far);

    let direction = ray.direction;
    let origin = ray.origin;

    let inverseDir = PointType.VectorType(repeating: 1.0) ./ direction;
    let unnormalizedT1 = (minBound - origin) .* inverseDir;
    let unnormalizedT2 = (maxBound - origin) .* inverseDir;
    let compareMask = unnormalizedT1 .> unnormalizedT2;
    let t1 = unnormalizedT1.replace(with: unnormalizedT2, where: compareMask)
    let t2 = unnormalizedT2.replace(with: unnormalizedT1, where: compareMask);
    tmin = PointType.VectorType.max(tmin, t1);
    tmax = PointType.VectorType.min(tmax, t2);
    if (tmin .> tmax).any {
      return nil;
    }
    return (tmin.maxElement(), tmax.minElement() + 0.01);
  }
}

public protocol Ray {
  associatedtype PointType: Point
  typealias ValueType = PointType.ValueType
  var direction: PointType.VectorType { get }
  var origin: PointType { get }
  
}

public protocol RenderContext {
  associatedtype ElementType
  associatedtype PointType where PointType.ValueType: IntConvertibleFloatingPoint,
                                 PointType == ElementType.PointType,
                                 PointType.AxisType == ElementType.PointType.MaskType.AxisType
  associatedtype CollisionType : Collision where CollisionType.ElementType == ElementType
  associatedtype RayType: Ray where RayType.PointType == PointType
  typealias VectorType = PointType.VectorType
  typealias ValueType = PointType.ValueType
  var bvhStack : BVHTree<ElementType>.IntersectionContext { get set }
}

let NUM_BUCKETS = 32
let MaxPrimitivesPerNode = 32
public enum HitMode {
  case AnyHit
  case NearestHit
}

public class BVHTree<Element: BVHElement>
where Element.PointType.ValueType : IntConvertibleFloatingPoint&Numeric,
      Element.PointType.AxisType == Element.PointType.MaskType.AxisType {
  typealias PointType = Element.PointType
  typealias VectorType = Element.PointType.VectorType
  public typealias ValueType = Element.PointType.ValueType
  typealias AxisType = Element.PointType.AxisType
  typealias BoundingBox = VectorTypes.BoundingBox<PointType>
  
  let root: Node
  // We rely on an owner existing for the primitives in the tree, so to
  // be safe we retain the array here
  let primitives: [Element]
  public class Node {
    init(_ bounds: BoundingBox, _ axis: AxisType, _ left: Node, _ right: Node) {
      isLeaf = false;
      self.axis = axis
      self.bounds = bounds
      self.children = (.passRetained(left), .passRetained(right))
      primitives = nil
    }
    init(_ bounds: BoundingBox, _ primitives: [Element.LeafElementStorage]) {
      self.axis = nil
      self.bounds = bounds
      isLeaf = true
      children = nil
      self.primitives = primitives
    }
    let isLeaf: Bool
    let bounds: BoundingBox
    let axis: AxisType?
    let children: (Unmanaged<Node>, Unmanaged<Node>)?
    let primitives: [Element.LeafElementStorage]?
  }
  private struct PrimitiveInfo {
    let element: Element
    let bounds: BoundingBox
    let centroid: PointType
  }
  private struct Bucket {
    var count = 0
    var leftInclusiveCount = 0
    var rightExclusiveCount = 0
    var leftInclusiveBounds = BoundingBox()
    var rightExclusiveBounds = BoundingBox()
    var leftInclusiveCentroidBounds = BoundingBox()
    var rightExclusiveCentroidBounds = BoundingBox()
    var bounds = BoundingBox()
    var centroidBounds = BoundingBox()
    var splitCost : ValueType = ValueType.zero
  }
  
  init(elements: [Element]) {
    self.primitives = elements
    var primitives = elements.map { element in
      PrimitiveInfo(element: element, bounds: element.bounds, centroid: element.bounds.centroid())
    }
    root = BVHTree.build(depth: 0, primitives: &primitives);
  }
  
  private static func makeLeaf<T: Sequence>(_ primitives: T) -> Node
  where T.Element == PrimitiveInfo
  {
    var bounds = BoundingBox()
    var list = [Element]()
//    var triangles = [Triangle]()
    for primitive in primitives {
      bounds = bounds.merge(other: primitive.bounds)
//      if let tri = primitive.element as? Triangle {
//        triangles.append(tri)
//        continue
//      }
      list.append(primitive.element)
    }
//    while triangles.count >= 4 {
//      list.append(.passRetained(TriangleBundle(triangles: triangles[0..<4].map({$0}))))
//      triangles.remove(at: 0)
//      triangles.remove(at: 0)
//      triangles.remove(at: 0)
//      triangles.remove(at: 0)
//    }
//    for t in triangles {
//      list.append(.passUnretained(t))
//    }
    return Node(bounds, Element.buildLeafElements(list))
  }
  private static func bucketForPrimitive(_ centroidBounds: BoundingBox,
                                         _ bucketCount: Int,
                                         _ axis: AxisType,
                                         _ primitive: PrimitiveInfo) -> Int {
    let relativePosition = centroidBounds.offsetRatio(point: primitive.centroid, axis: axis);
    return min((ValueType(NUM_BUCKETS) * relativePosition).asInteger, bucketCount - 1);
  }
  private static func build(depth: Int, primitives: inout [PrimitiveInfo]) -> Node {
    assert(!primitives.isEmpty)
    var bounds = BoundingBox()
    for primitive in primitives {
      bounds = bounds.merge(other: primitive.bounds);
    }

    if primitives.count <= 4 || bounds.surfaceArea() == 0{
      return makeLeaf(primitives)
    }
    var centroidBounds = BoundingBox()
    for primitive in primitives {
      centroidBounds = centroidBounds.merge(point: primitive.centroid)
    }
    let maxAxis = centroidBounds.maxAxis();
    if (centroidBounds.maxBound[maxAxis] == centroidBounds.minBound[maxAxis]) {
      return makeLeaf(primitives)
    }
    var buckets = [Bucket](repeating: Bucket(), count: NUM_BUCKETS)
    
    // First pass, accrue the per bucket primitive information
    for primitive in primitives {
      let b = bucketForPrimitive(centroidBounds, NUM_BUCKETS, maxAxis, primitive);
      assert(b < NUM_BUCKETS);
      buckets[b].count += 1;
      buckets[b].bounds = buckets[b].bounds.merge(other: primitive.bounds);
      buckets[b].centroidBounds = buckets[b].centroidBounds.merge(point: primitive.centroid);
    }
    
    // Second pass, compute left properties
    do {
      var cummulativeBounds = BoundingBox();
      var cummulativeCentroidBounds = BoundingBox();
      var cummulativeCount = 0;
      for i in 0..<NUM_BUCKETS {
        buckets[i].leftInclusiveCount = cummulativeCount + buckets[i].count;
        cummulativeCount = buckets[i].leftInclusiveCount;
        buckets[i].leftInclusiveBounds = cummulativeBounds.merge(other: buckets[i].bounds);
        cummulativeBounds = buckets[i].leftInclusiveBounds;
        buckets[i].leftInclusiveCentroidBounds = cummulativeCentroidBounds.merge(other: buckets[i].leftInclusiveCentroidBounds);
        cummulativeCentroidBounds = buckets[i].leftInclusiveCentroidBounds;
      }
    }
    // Third pass: compute right properties
    do {
      var priorBounds = BoundingBox();
      var priorCentroidBounds = BoundingBox();
      var cummulativeCount = 0;
      for j in 0..<NUM_BUCKETS {
        let i = NUM_BUCKETS - j - 1;
        buckets[i].rightExclusiveBounds = priorBounds;
        priorBounds = priorBounds.merge(other: buckets[i].bounds);
        buckets[i].rightExclusiveCentroidBounds = priorCentroidBounds;
        priorCentroidBounds = priorCentroidBounds.merge(other: buckets[i].rightExclusiveCentroidBounds);
        buckets[i].rightExclusiveCount = cummulativeCount;
        cummulativeCount += buckets[i].count;
      }
    }
    // Fourth pass: compute split costs
    do {
      let leafSurface = bounds.surfaceArea();
      assert(leafSurface>0)
      for i in 0..<buckets.count {
        assert(buckets[i].leftInclusiveCount + buckets[i].rightExclusiveCount == primitives.count);
        let leftCost : ValueType = buckets[i].leftInclusiveBounds.surfaceArea() / leafSurface * ValueType(buckets[i].leftInclusiveCount);
        let rightCost : ValueType =
          buckets[i].rightExclusiveBounds.surfaceArea() / leafSurface * ValueType(buckets[i].rightExclusiveCount);
        buckets[i].splitCost = 1.0 + leftCost * 2.0 + rightCost * 2.0;
      }
    }
    
    // Find the cheapest split
    var minimumSplitCost = buckets[0].splitCost
    var minimumSplitIndex = 0
    for i in 1..<buckets.count {
      if buckets[i].splitCost < minimumSplitCost {
        minimumSplitCost = buckets[i].splitCost
        minimumSplitIndex = i
      }
      
    }
    
    let leafCost = ValueType(primitives.count)
    if leafCost < minimumSplitCost && primitives.count < MaxPrimitivesPerNode {
      return makeLeaf(primitives)
    }
    let firstOverPivot = primitives.partition { info in
      let dest = bucketForPrimitive(centroidBounds, NUM_BUCKETS, maxAxis, info)
      return dest > minimumSplitIndex
    }
    var leftPrimitives = Array(primitives[primitives.startIndex..<firstOverPivot])
    var rightPrimitives = Array(primitives[firstOverPivot...])
    let left = build(depth: depth + 1, primitives: &leftPrimitives)
    let right = build(depth: depth + 1, primitives: &rightPrimitives)
    return Node(bounds, maxAxis, left, right)
  }

  func intersect<RenderContextType: RenderContext>(
    context: RenderContextType,
    ray: RenderContextType.RayType,
    hitMode: HitMode,
    min: ValueType,
    max: ValueType
  ) -> RenderContextType.CollisionType?
  where RenderContextType.ElementType == Element {
    return intersect(context: context, ray, hitMode, min, max);
  }
  public class Stack<T> {
    var currentIndex = 0
    var buffer = [T]()
    @inline(__always) func push(_ t: T) {
      if buffer.count <= currentIndex {
        buffer.append(t)
      } else {
        buffer[currentIndex] = t
      }
      currentIndex += 1
    }
    @inline(__always) func pop() -> T? {
      if currentIndex == 0 {
        return nil
      }
      currentIndex -= 1
      return buffer[currentIndex]
    }
    func reset() {
      currentIndex = 0
    }
  }
  
  public typealias IntersectionContext = Stack<(Unmanaged<Node>, min: ValueType, max: ValueType)>

  private func intersect<RenderContextType: RenderContext> (
    context: RenderContextType,
    _ ray: RenderContextType.RayType,
    _ hitMode: HitMode,
    _ parentMin: ValueType,
    _ parentMax: ValueType
  ) -> RenderContextType.CollisionType?
  where RenderContextType.ElementType == Element {
    context.bvhStack.push((.passUnretained(root), parentMin, parentMax));
    var result : RenderContextType.CollisionType? = nil;
    var nearest = parentMax;
    var primitiveCount = 0;
    var nodeCount = 0;
  traversal_loop: while let node = context.bvhStack.pop() {
      let (value, nodeMin, nodeMax) = node;
      nodeCount += 1;
      let dirIsNegative = ray.direction .< VectorType(repeating: 0)
      if nodeMin > nearest {
        continue;
      }

      let farIntersect = min(nearest, nodeMax);
      guard !value._withUnsafeGuaranteedRef({$0.isLeaf}) else {
        guard let collisionBounds = value._withUnsafeGuaranteedRef({$0.bounds}).intersect(ray, min: nodeMin - 0.01, max: nearest) else {
          continue traversal_loop
        }

        primitiveCount += value._withUnsafeGuaranteedRef { $0.primitives!.count }
        for i in 0..<value._withUnsafeGuaranteedRef({ $0.primitives! }).count { // retain
          let primitive = value._withUnsafeGuaranteedRef({ $0.primitives![i] })
          guard let collision = primitive.intersect(context: context, ray: ray, hitMode: hitMode, min: max(nodeMin, collisionBounds.min), max: min(nearest, collisionBounds.max)) else {
            continue
          }
          if collision.distance < nearest {
            nearest = collision.distance;
            result = collision
            if hitMode == .AnyHit {
              break traversal_loop
            }
          }
        }
        continue traversal_loop
      }

      guard let boundsCollision = value._withUnsafeGuaranteedRef({$0.bounds}).intersect(ray, min: nodeMin, max: farIntersect) else {
        continue
      }
      let (childMin, childMax) = boundsCollision;
      if dirIsNegative[value._withUnsafeGuaranteedRef({$0.axis!})] {
        context.bvhStack.push((value._withUnsafeGuaranteedRef({$0.children!.1}), childMin, childMax));
        context.bvhStack.push((value._withUnsafeGuaranteedRef({$0.children!.0}), childMin, childMax));
      } else {
        context.bvhStack.push((value._withUnsafeGuaranteedRef({$0.children!.0}), childMin, childMax));
        context.bvhStack.push((value._withUnsafeGuaranteedRef({$0.children!.1}), childMin, childMax));
      }
    }

    if let c = result {
      result = RenderContextType.CollisionType(distance: c.distance,
                                               uv: c.uv,
                                               intersectionCount: primitiveCount,
                                               nodeCount: nodeCount,
                                               target: c.target)
    }
    return result;
  }
}

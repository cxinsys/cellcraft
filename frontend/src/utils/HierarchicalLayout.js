/**
 * 커스텀 계층적 레이아웃 알고리즘
 * Dagre 없이 토폴로지컬 정렬을 기반으로 DAG 레이아웃을 구성
 */

export class HierarchicalLayout {
  constructor(nodes, edges, options = {}) {
    this.nodes = nodes;
    this.edges = edges;
    this.options = {
      direction: options.direction || 'TB', // TB, LR, BT, RL
      rankSep: options.rankSep || 180,      // 레벨 간 거리
      nodeSep: options.nodeSep || 250,      // 노드 간 거리
      marginX: options.marginX || 50,
      marginY: options.marginY || 50,
      ...options
    };
    
    this.levels = {};
    this.positions = {};
    this.levelNodes = {};
  }

  /**
   * 토폴로지컬 정렬을 기반으로 레벨 할당
   */
  assignLevels() {
    const inDegree = {};
    const adjList = {};
    const visited = new Set();
    
    // 초기화
    this.nodes.forEach(node => {
      inDegree[node.id] = 0;
      adjList[node.id] = [];
    });
    
    // 인접 리스트와 진입 차수 구성
    this.edges.forEach(edge => {
      adjList[edge.source].push(edge.target);
      inDegree[edge.target]++;
    });
    
    // BFS로 레벨 할당
    const queue = [];
    
    // 진입 차수가 0인 노드들로 시작
    this.nodes.forEach(node => {
      if (inDegree[node.id] === 0) {
        queue.push(node.id);
        this.levels[node.id] = 0;
      }
    });
    
    while (queue.length > 0) {
      const current = queue.shift();
      visited.add(current);
      
      // 인접한 노드들 처리
      adjList[current].forEach(neighbor => {
        inDegree[neighbor]--;
        
        if (inDegree[neighbor] === 0) {
          // 최대 부모 레벨 + 1로 설정
          const parentLevels = this.edges
            .filter(edge => edge.target === neighbor)
            .map(edge => this.levels[edge.source])
            .filter(level => level !== undefined);
          
          this.levels[neighbor] = parentLevels.length > 0 ? 
            Math.max(...parentLevels) + 1 : 0;
          
          queue.push(neighbor);
        }
      });
    }
    
    // 순환 참조가 있는 경우 처리
    this.nodes.forEach(node => {
      if (!visited.has(node.id)) {
        this.levels[node.id] = 0;
      }
    });
    
    return this.levels;
  }

  /**
   * 방향에 따른 위치 계산
   */
  calculatePositions(direction = null) {
    const dir = direction || this.options.direction;
    this.assignLevels();
    
    // 레벨별 노드 그룹화
    this.levelNodes = {};
    Object.entries(this.levels).forEach(([nodeId, level]) => {
      if (!this.levelNodes[level]) {
        this.levelNodes[level] = [];
      }
      this.levelNodes[level].push(nodeId);
    });
    
    // 각 레벨에서 엣지 교차 최소화
    this.minimizeCrossings();
    
    // 위치 계산
    Object.entries(this.levels).forEach(([nodeId, level]) => {
      const levelNodes = this.levelNodes[level];
      const index = levelNodes.indexOf(nodeId);
      const count = levelNodes.length;
      
      // 동적 간격 조정
      const dynamicSep = this.calculateDynamicSpacing(count);
      
      // 중앙 정렬을 위한 오프셋 계산
      const offset = (index - (count - 1) / 2) * dynamicSep;
      
      // 방향에 따른 좌표 설정
      let x, y;
      switch(dir) {
        case 'TB': // Top to Bottom
          x = offset;
          y = level * this.options.rankSep;
          break;
        case 'BT': // Bottom to Top
          x = offset;
          y = -level * this.options.rankSep;
          break;
        case 'LR': // Left to Right
          x = level * this.options.rankSep;
          y = offset;
          break;
        case 'RL': // Right to Left
          x = -level * this.options.rankSep;
          y = offset;
          break;
        default:
          x = offset;
          y = level * this.options.rankSep;
      }
      
      this.positions[nodeId] = { 
        x, 
        y, 
        level,
        index: index,
        levelNodeCount: count
      };
    });
    
    return this.positions;
  }

  /**
   * 동적 간격 조정
   */
  calculateDynamicSpacing(nodeCount) {
    const baseSpacing = this.options.nodeSep;
    const minSpacing = 150;
    const maxSpacing = 400;
    
    // 노드 수가 많을 때는 간격을 줄이고, 적을 때는 늘림
    let spacing = baseSpacing;
    
    if (nodeCount > 4) {
      spacing = Math.max(minSpacing, baseSpacing * (4 / nodeCount));
    } else if (nodeCount < 3) {
      spacing = Math.min(maxSpacing, baseSpacing * 1.2);
    }
    
    return spacing;
  }

  /**
   * 바리센터 휴리스틱을 사용한 엣지 교차 최소화
   */
  minimizeCrossings() {
    const maxIterations = 10;
    let improved = true;
    let iteration = 0;
    
    while (improved && iteration < maxIterations) {
      improved = false;
      iteration++;
      
      // 각 레벨에서 노드 재정렬
      Object.keys(this.levelNodes)
        .sort((a, b) => parseInt(a) - parseInt(b))
        .forEach(level => {
          if (parseInt(level) > 0) {
            const oldOrder = [...this.levelNodes[level]];
            this.reorderLevel(parseInt(level));
            
            // 순서가 변경되었는지 확인
            if (!this.arraysEqual(oldOrder, this.levelNodes[level])) {
              improved = true;
            }
          }
        });
    }
  }

  /**
   * 레벨 내 노드 재정렬 (바리센터 휴리스틱)
   */
  reorderLevel(currentLevel) {
    const currentLevelNodes = this.levelNodes[currentLevel];
    const previousLevel = currentLevel - 1;
    
    if (!this.levelNodes[previousLevel]) return;
    
    const previousLevelNodes = this.levelNodes[previousLevel];
    
    // 각 노드의 바리센터 값 계산
    const barycenters = {};
    
    currentLevelNodes.forEach(nodeId => {
      const parents = this.edges
        .filter(edge => edge.target === nodeId && 
                       previousLevelNodes.includes(edge.source))
        .map(edge => previousLevelNodes.indexOf(edge.source));
      
      if (parents.length > 0) {
        barycenters[nodeId] = parents.reduce((a, b) => a + b, 0) / parents.length;
      } else {
        // 부모가 없는 경우 현재 위치 유지
        barycenters[nodeId] = currentLevelNodes.indexOf(nodeId);
      }
    });
    
    // 바리센터 값으로 정렬
    currentLevelNodes.sort((a, b) => {
      const aCenter = barycenters[a] !== undefined ? barycenters[a] : 
                     currentLevelNodes.indexOf(a);
      const bCenter = barycenters[b] !== undefined ? barycenters[b] : 
                     currentLevelNodes.indexOf(b);
      
      return aCenter - bCenter;
    });
  }

  /**
   * 배열 동일성 검사
   */
  arraysEqual(arr1, arr2) {
    if (arr1.length !== arr2.length) return false;
    return arr1.every((val, index) => val === arr2[index]);
  }

  /**
   * 특수 케이스 처리 (TENET 파이프라인 최적화)
   */
  handleSpecialCases() {
    // GRN_FDR과 GRN_NumLinks가 같은 레벨에 있는 경우
    const level2Nodes = this.levelNodes[2];
    if (level2Nodes && level2Nodes.length >= 2) {
      const grnFdrIndex = level2Nodes.findIndex(id => 
        id.includes('GRN') && id.includes('FDR'));
      const grnNumIndex = level2Nodes.findIndex(id => 
        id.includes('GRN') && id.includes('NumLinks'));
      
      if (grnFdrIndex !== -1 && grnNumIndex !== -1) {
        // GRN_FDR을 왼쪽에, GRN_NumLinks를 오른쪽에 배치
        if (grnFdrIndex > grnNumIndex) {
          [level2Nodes[grnFdrIndex], level2Nodes[grnNumIndex]] = 
          [level2Nodes[grnNumIndex], level2Nodes[grnFdrIndex]];
        }
      }
    }
  }

  /**
   * 레이아웃 방향 변경
   */
  changeDirection(newDirection) {
    this.options.direction = newDirection;
    return this.calculatePositions(newDirection);
  }

  /**
   * 컨테이너 크기에 맞춰 스케일링
   */
  scaleToFit(containerWidth, containerHeight, padding = 50) {
    if (Object.keys(this.positions).length === 0) {
      this.calculatePositions();
    }
    
    // 현재 레이아웃의 경계 계산
    const xValues = Object.values(this.positions).map(pos => pos.x);
    const yValues = Object.values(this.positions).map(pos => pos.y);
    
    const minX = Math.min(...xValues);
    const maxX = Math.max(...xValues);
    const minY = Math.min(...yValues);
    const maxY = Math.max(...yValues);
    
    const currentWidth = maxX - minX;
    const currentHeight = maxY - minY;
    
    // 사용 가능한 공간 계산 (패딩 제외)
    const availableWidth = containerWidth - (padding * 2);
    const availableHeight = containerHeight - (padding * 2);
    
    // 스케일 팩터 계산
    const scaleX = currentWidth > 0 ? availableWidth / currentWidth : 1;
    const scaleY = currentHeight > 0 ? availableHeight / currentHeight : 1;
    const scale = Math.min(scaleX, scaleY, 1); // 최대 1배까지만 확대
    
    // 위치 조정 (중앙 정렬)
    const centerX = (minX + maxX) / 2;
    const centerY = (minY + maxY) / 2;
    
    const scaledPositions = {};
    Object.entries(this.positions).forEach(([nodeId, pos]) => {
      scaledPositions[nodeId] = {
        ...pos,
        x: (pos.x - centerX) * scale,
        y: (pos.y - centerY) * scale
      };
    });
    
    return {
      positions: scaledPositions,
      scale: scale,
      bounds: {
        width: currentWidth * scale,
        height: currentHeight * scale,
        minX: (minX - centerX) * scale,
        maxX: (maxX - centerX) * scale,
        minY: (minY - centerY) * scale,
        maxY: (maxY - centerY) * scale
      }
    };
  }

  /**
   * Critical Path 찾기 (최장 경로)
   */
  findCriticalPath() {
    const distances = {};
    const parents = {};
    
    // 토폴로지컬 정렬된 순서로 처리
    const sorted = this.topologicalSort();
    
    // 초기화
    sorted.forEach(nodeId => {
      distances[nodeId] = 0;
      parents[nodeId] = null;
    });
    
    // 최장 경로 계산
    sorted.forEach(nodeId => {
      this.edges
        .filter(edge => edge.source === nodeId)
        .forEach(edge => {
          const target = edge.target;
          if (distances[target] < distances[nodeId] + 1) {
            distances[target] = distances[nodeId] + 1;
            parents[target] = nodeId;
          }
        });
    });
    
    // 최장 경로 역추적
    let endNode = Object.keys(distances).reduce((a, b) => 
      distances[a] > distances[b] ? a : b
    );
    
    const path = [];
    while (endNode !== null) {
      path.unshift(endNode);
      endNode = parents[endNode];
    }
    
    return path;
  }

  /**
   * 토폴로지컬 정렬
   */
  topologicalSort() {
    const inDegree = {};
    const queue = [];
    const result = [];
    
    // 진입 차수 계산
    this.nodes.forEach(node => inDegree[node.id] = 0);
    this.edges.forEach(edge => inDegree[edge.target]++);
    
    // 시작 노드들 큐에 추가
    Object.keys(inDegree).forEach(nodeId => {
      if (inDegree[nodeId] === 0) {
        queue.push(nodeId);
      }
    });
    
    // BFS 실행
    while (queue.length > 0) {
      const current = queue.shift();
      result.push(current);
      
      this.edges
        .filter(edge => edge.source === current)
        .forEach(edge => {
          inDegree[edge.target]--;
          if (inDegree[edge.target] === 0) {
            queue.push(edge.target);
          }
        });
    }
    
    return result;
  }

  /**
   * 레이아웃 정보 반환
   */
  getLayoutInfo() {
    return {
      nodes: this.nodes,
      edges: this.edges,
      positions: this.positions,
      levels: this.levels,
      levelNodes: this.levelNodes,
      options: this.options,
      criticalPath: this.findCriticalPath()
    };
  }
}

export default HierarchicalLayout;
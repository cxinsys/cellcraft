/**
 * MSW (Mock Service Worker) API 핸들러
 *
 * 테스트 환경에서 사용할 API 모킹 정의
 * 실제 백엔드 호출 대신 모의 응답을 반환합니다.
 */

import { http, HttpResponse } from 'msw';

// 기본 API URL (환경에 따라 변경)
const API_BASE = 'http://localhost:8000/api';

export const handlers = [
  /**
   * 로그인 API
   * POST /routes/auth/login/access-token
   */
  http.post(`${API_BASE}/routes/auth/login/access-token`, async ({ request }) => {
    const formData = await request.formData();
    const username = formData.get('username');
    const password = formData.get('password');

    // 올바른 자격증명
    if (username === 'testuser' && password === 'testpassword') {
      return HttpResponse.json({
        access_token: 'mock-jwt-token',
        token_type: 'bearer'
      });
    }

    // 잘못된 자격증명
    return HttpResponse.json(
      { detail: 'Incorrect username or password' },
      { status: 401 }
    );
  }),

  /**
   * 사용자 정보 조회
   * GET /routes/user/me
   */
  http.get(`${API_BASE}/routes/user/me`, ({ request }) => {
    const authHeader = request.headers.get('Authorization');

    // 인증되지 않은 요청
    if (!authHeader || !authHeader.includes('Bearer')) {
      return HttpResponse.json(
        { detail: 'Not authenticated' },
        { status: 401 }
      );
    }

    // 인증된 사용자 정보 반환
    return HttpResponse.json({
      id: 1,
      username: 'testuser',
      email: 'test@example.com',
      is_superuser: false
    });
  }),

  /**
   * 워크플로우 목록 조회
   * GET /routes/workflow
   */
  http.get(`${API_BASE}/routes/workflow`, () => {
    return HttpResponse.json([
      {
        id: 1,
        name: 'Test Workflow 1',
        created_at: '2024-01-01T00:00:00Z',
        updated_at: '2024-01-01T00:00:00Z'
      },
      {
        id: 2,
        name: 'Test Workflow 2',
        created_at: '2024-01-02T00:00:00Z',
        updated_at: '2024-01-02T00:00:00Z'
      }
    ]);
  }),

  /**
   * 워크플로우 생성
   * POST /routes/workflow
   */
  http.post(`${API_BASE}/routes/workflow`, async ({ request }) => {
    const body = await request.json();

    return HttpResponse.json({
      id: 3,
      name: body.name || 'New Workflow',
      created_at: new Date().toISOString(),
      updated_at: new Date().toISOString()
    }, { status: 201 });
  }),

  /**
   * 파일 업로드
   * POST /api/upload
   */
  http.post(`${API_BASE}/upload`, async ({ request }) => {
    const formData = await request.formData();
    const file = formData.get('file');

    if (!file) {
      return HttpResponse.json(
        { detail: 'No file provided' },
        { status: 400 }
      );
    }

    return HttpResponse.json({
      filename: file.name,
      path: `/uploads/${file.name}`,
      size: file.size
    });
  }),

  /**
   * 플러그인 목록 조회
   * GET /routes/plugin
   */
  http.get(`${API_BASE}/routes/plugin`, () => {
    return HttpResponse.json([
      {
        plugin_name: 'TENET',
        display_name: 'TENET',
        description: 'Gene Regulatory Network analysis',
        version: '1.0.0'
      },
      {
        plugin_name: 'GENIE3',
        display_name: 'GENIE3',
        description: 'Gene network inference',
        version: '1.0.0'
      }
    ]);
  }),

  /**
   * 태스크 제출
   * POST /api/task/submit
   */
  http.post(`${API_BASE}/task/submit`, async ({ request }) => {
    const body = await request.json();

    return HttpResponse.json({
      task_id: 'mock-task-id-123',
      status: 'submitted',
      workflow_id: body.workflow_id || null
    });
  }),

  /**
   * 태스크 상태 조회
   * GET /api/task/:taskId/status
   */
  http.get(`${API_BASE}/task/:taskId/status`, ({ params }) => {
    return HttpResponse.json({
      task_id: params.taskId,
      status: 'completed',
      progress: 100,
      result: {
        output_files: ['result1.csv', 'result2.png']
      }
    });
  }),
];

/**
 * 테스트에서 사용하는 방법:
 *
 * 1. Vitest에서 사용:
 *    import { setupServer } from 'msw/node';
 *    import { handlers } from './tests/mocks/handlers';
 *
 *    const server = setupServer(...handlers);
 *    beforeAll(() => server.listen());
 *    afterEach(() => server.resetHandlers());
 *    afterAll(() => server.close());
 *
 * 2. 개별 테스트에서 핸들러 오버라이드:
 *    server.use(
 *      http.get('/api/routes/workflow', () => {
 *        return HttpResponse.json([]);
 *      })
 *    );
 */

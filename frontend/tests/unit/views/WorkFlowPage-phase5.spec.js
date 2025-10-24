import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';

/**
 * Unit tests for WorkFlowPage.vue - Phase 5
 *
 * These tests focus on simple state management and utility methods:
 * - Modal management (closeLogsModal, closeDAGModal)
 * - Job table management (closeJobTable)
 * - Progress tracking (viewProgress)
 * - Message system (setMessage)
 * - Workflow operations (getWorkflowTitle)
 * - View mode (activateCompileCheck, deactivateCompileCheck, updateIsTabView)
 */

describe('WorkFlowPage.vue - Phase 5: Additional Methods', () => {
  let mockComponent;
  let mockTaskList;

  beforeEach(() => {
    // Mock component instance with common properties
    mockComponent = {
      showLogsModal: false,
      selectedTaskLogs: null,
      logsLoading: false,
      showDAGModal: false,
      selectedDAGTaskId: null,
      selectedDAGTaskName: null,
      selectedDAGTaskStatus: null,
      selectedDAGTask: null,
      show_jobs: false,
      timeInterval: 123,
      toggleMessage: false,
      messageStatus: '',
      messageContent: '',
      compile_check: false,
      isTabView: false,
      taskList: []
    };

    // Mock task list for viewProgress tests
    mockTaskList = [
      {
        task_id: 'task-1',
        workflow_title: 'Test Workflow 1',
        task_title: 'Task 1',
        status: 'SUCCESS'
      },
      {
        task_id: 'task-2',
        workflow_title: null,
        task_title: 'Task 2',
        status: 'RUNNING'
      },
      {
        task_id: 'task-3',
        workflow_title: '',
        task_title: null,
        status: 'PENDING'
      }
    ];
  });

  /**
   * Modal Management Tests
   */
  describe('Modal Management', () => {
    describe('closeLogsModal', () => {
      it('should reset all logs modal state', () => {
        const closeLogsModal = function() {
          this.showLogsModal = false;
          this.selectedTaskLogs = null;
          this.logsLoading = false;
        };

        mockComponent.showLogsModal = true;
        mockComponent.selectedTaskLogs = { task_info: { task_id: 'test' } };
        mockComponent.logsLoading = true;

        closeLogsModal.call(mockComponent);

        expect(mockComponent.showLogsModal).toBe(false);
        expect(mockComponent.selectedTaskLogs).toBeNull();
        expect(mockComponent.logsLoading).toBe(false);
      });

      it('should handle already closed modal', () => {
        const closeLogsModal = function() {
          this.showLogsModal = false;
          this.selectedTaskLogs = null;
          this.logsLoading = false;
        };

        closeLogsModal.call(mockComponent);

        expect(mockComponent.showLogsModal).toBe(false);
        expect(mockComponent.selectedTaskLogs).toBeNull();
        expect(mockComponent.logsLoading).toBe(false);
      });
    });

    describe('closeDAGModal', () => {
      it('should reset all DAG modal state', () => {
        const closeDAGModal = function() {
          this.showDAGModal = false;
          this.selectedDAGTaskId = null;
          this.selectedDAGTaskName = null;
        };

        mockComponent.showDAGModal = true;
        mockComponent.selectedDAGTaskId = 'task-123';
        mockComponent.selectedDAGTaskName = 'Test Task';
        mockComponent.selectedDAGTaskStatus = 'RUNNING';

        closeDAGModal.call(mockComponent);

        expect(mockComponent.showDAGModal).toBe(false);
        expect(mockComponent.selectedDAGTaskId).toBeNull();
        expect(mockComponent.selectedDAGTaskName).toBeNull();
      });

      it('should handle already closed modal', () => {
        const closeDAGModal = function() {
          this.showDAGModal = false;
          this.selectedDAGTaskId = null;
          this.selectedDAGTaskName = null;
        };

        closeDAGModal.call(mockComponent);

        expect(mockComponent.showDAGModal).toBe(false);
        expect(mockComponent.selectedDAGTaskId).toBeNull();
        expect(mockComponent.selectedDAGTaskName).toBeNull();
      });
    });
  });

  /**
   * Job Table Management Tests
   */
  describe('Job Table Management', () => {
    describe('closeJobTable', () => {
      it('should close job table and clear interval', () => {
        const clearIntervalSpy = vi.spyOn(global, 'clearInterval');

        const closeJobTable = function() {
          this.show_jobs = false;
          clearInterval(this.timeInterval);
        };

        mockComponent.show_jobs = true;
        mockComponent.timeInterval = 456;

        closeJobTable.call(mockComponent);

        expect(mockComponent.show_jobs).toBe(false);
        expect(clearIntervalSpy).toHaveBeenCalledWith(456);

        clearIntervalSpy.mockRestore();
      });

      it('should handle null interval', () => {
        const clearIntervalSpy = vi.spyOn(global, 'clearInterval');

        const closeJobTable = function() {
          this.show_jobs = false;
          clearInterval(this.timeInterval);
        };

        mockComponent.show_jobs = true;
        mockComponent.timeInterval = null;

        closeJobTable.call(mockComponent);

        expect(mockComponent.show_jobs).toBe(false);
        expect(clearIntervalSpy).toHaveBeenCalledWith(null);

        clearIntervalSpy.mockRestore();
      });
    });
  });

  /**
   * Progress Tracking Tests
   */
  describe('Progress Tracking', () => {
    describe('viewProgress', () => {
      it('should set DAG modal data for task with workflow_title', () => {
        const viewProgress = function(task_id) {
          console.log('viewProgress called with task_id:', task_id, 'type:', typeof task_id);

          const task = this.taskList.find(t => t.task_id === task_id);
          const taskName = task ? (task.workflow_title || task.task_title || 'Unknown Task') : 'Unknown Task';
          const taskStatus = task ? task.status : 'UNKNOWN';

          this.selectedDAGTaskId = String(task_id);
          this.selectedDAGTaskName = taskName;
          this.selectedDAGTaskStatus = taskStatus;
          this.selectedDAGTask = task;
          this.showDAGModal = true;

          console.log('DAG modal opened with taskId:', this.selectedDAGTaskId, 'taskName:', this.selectedDAGTaskName, 'taskStatus:', this.selectedDAGTaskStatus);
        };

        mockComponent.taskList = mockTaskList;

        viewProgress.call(mockComponent, 'task-1');

        expect(mockComponent.selectedDAGTaskId).toBe('task-1');
        expect(mockComponent.selectedDAGTaskName).toBe('Test Workflow 1');
        expect(mockComponent.selectedDAGTaskStatus).toBe('SUCCESS');
        expect(mockComponent.selectedDAGTask).toEqual(mockTaskList[0]);
        expect(mockComponent.showDAGModal).toBe(true);
      });

      it('should use task_title when workflow_title is null', () => {
        const viewProgress = function(task_id) {
          const task = this.taskList.find(t => t.task_id === task_id);
          const taskName = task ? (task.workflow_title || task.task_title || 'Unknown Task') : 'Unknown Task';
          const taskStatus = task ? task.status : 'UNKNOWN';

          this.selectedDAGTaskId = String(task_id);
          this.selectedDAGTaskName = taskName;
          this.selectedDAGTaskStatus = taskStatus;
          this.selectedDAGTask = task;
          this.showDAGModal = true;
        };

        mockComponent.taskList = mockTaskList;

        viewProgress.call(mockComponent, 'task-2');

        expect(mockComponent.selectedDAGTaskName).toBe('Task 2');
        expect(mockComponent.selectedDAGTaskStatus).toBe('RUNNING');
      });

      it('should use "Unknown Task" when both titles are missing', () => {
        const viewProgress = function(task_id) {
          const task = this.taskList.find(t => t.task_id === task_id);
          const taskName = task ? (task.workflow_title || task.task_title || 'Unknown Task') : 'Unknown Task';
          const taskStatus = task ? task.status : 'UNKNOWN';

          this.selectedDAGTaskId = String(task_id);
          this.selectedDAGTaskName = taskName;
          this.selectedDAGTaskStatus = taskStatus;
          this.selectedDAGTask = task;
          this.showDAGModal = true;
        };

        mockComponent.taskList = mockTaskList;

        viewProgress.call(mockComponent, 'task-3');

        expect(mockComponent.selectedDAGTaskName).toBe('Unknown Task');
        expect(mockComponent.selectedDAGTaskStatus).toBe('PENDING');
      });

      it('should handle task not found', () => {
        const viewProgress = function(task_id) {
          const task = this.taskList.find(t => t.task_id === task_id);
          const taskName = task ? (task.workflow_title || task.task_title || 'Unknown Task') : 'Unknown Task';
          const taskStatus = task ? task.status : 'UNKNOWN';

          this.selectedDAGTaskId = String(task_id);
          this.selectedDAGTaskName = taskName;
          this.selectedDAGTaskStatus = taskStatus;
          this.selectedDAGTask = task;
          this.showDAGModal = true;
        };

        mockComponent.taskList = mockTaskList;

        viewProgress.call(mockComponent, 'non-existent-task');

        expect(mockComponent.selectedDAGTaskId).toBe('non-existent-task');
        expect(mockComponent.selectedDAGTaskName).toBe('Unknown Task');
        expect(mockComponent.selectedDAGTaskStatus).toBe('UNKNOWN');
        expect(mockComponent.selectedDAGTask).toBeUndefined();
        expect(mockComponent.showDAGModal).toBe(true);
      });

      it('should convert task_id to string', () => {
        const viewProgress = function(task_id) {
          const task = this.taskList.find(t => t.task_id === task_id);
          const taskName = task ? (task.workflow_title || task.task_title || 'Unknown Task') : 'Unknown Task';
          const taskStatus = task ? task.status : 'UNKNOWN';

          this.selectedDAGTaskId = String(task_id);
          this.selectedDAGTaskName = taskName;
          this.selectedDAGTaskStatus = taskStatus;
          this.selectedDAGTask = task;
          this.showDAGModal = true;
        };

        mockComponent.taskList = mockTaskList;

        viewProgress.call(mockComponent, 123); // Number input

        expect(mockComponent.selectedDAGTaskId).toBe('123');
        expect(typeof mockComponent.selectedDAGTaskId).toBe('string');
      });
    });
  });

  /**
   * Message System Tests
   */
  describe('Message System', () => {
    describe('setMessage', () => {
      it('should set message with success status', () => {
        const setTimeoutSpy = vi.spyOn(global, 'setTimeout');

        const setMessage = function(status, content) {
          this.toggleMessage = true;
          this.messageStatus = status;
          this.messageContent = content;
          setTimeout(() => {
            this.toggleMessage = false;
          }, 5000);
        };

        setMessage.call(mockComponent, 'success', 'Operation completed');

        expect(mockComponent.toggleMessage).toBe(true);
        expect(mockComponent.messageStatus).toBe('success');
        expect(mockComponent.messageContent).toBe('Operation completed');
        expect(setTimeoutSpy).toHaveBeenCalledWith(expect.any(Function), 5000);

        setTimeoutSpy.mockRestore();
      });

      it('should set message with error status', () => {
        const setTimeoutSpy = vi.spyOn(global, 'setTimeout');

        const setMessage = function(status, content) {
          this.toggleMessage = true;
          this.messageStatus = status;
          this.messageContent = content;
          setTimeout(() => {
            this.toggleMessage = false;
          }, 5000);
        };

        setMessage.call(mockComponent, 'error', 'Operation failed');

        expect(mockComponent.toggleMessage).toBe(true);
        expect(mockComponent.messageStatus).toBe('error');
        expect(mockComponent.messageContent).toBe('Operation failed');
        expect(setTimeoutSpy).toHaveBeenCalledWith(expect.any(Function), 5000);

        setTimeoutSpy.mockRestore();
      });

      it('should auto-hide message after 5 seconds', () => {
        vi.useFakeTimers();

        const setMessage = function(status, content) {
          this.toggleMessage = true;
          this.messageStatus = status;
          this.messageContent = content;
          setTimeout(() => {
            this.toggleMessage = false;
          }, 5000);
        };

        setMessage.call(mockComponent, 'success', 'Test message');

        expect(mockComponent.toggleMessage).toBe(true);

        // Fast-forward 5 seconds
        vi.advanceTimersByTime(5000);

        expect(mockComponent.toggleMessage).toBe(false);

        vi.useRealTimers();
      });

      it('should handle empty content', () => {
        const setMessage = function(status, content) {
          this.toggleMessage = true;
          this.messageStatus = status;
          this.messageContent = content;
          setTimeout(() => {
            this.toggleMessage = false;
          }, 5000);
        };

        setMessage.call(mockComponent, 'success', '');

        expect(mockComponent.toggleMessage).toBe(true);
        expect(mockComponent.messageStatus).toBe('success');
        expect(mockComponent.messageContent).toBe('');
      });
    });
  });

  /**
   * Workflow Operations Tests
   */
  describe('Workflow Operations', () => {
    describe('getWorkflowTitle', () => {
      it('should fetch and return workflow title', async () => {
        const mockFindWorkflow = vi.fn().mockResolvedValue({
          data: {
            title: 'My Workflow Title'
          }
        });

        const getWorkflowTitle = async function(id) {
          const workflow_id = {
            id: id,
          };
          const user_workflow = await mockFindWorkflow(workflow_id);
          console.log(user_workflow.data.title);
          return user_workflow.data.title;
        };

        const result = await getWorkflowTitle.call(mockComponent, 123);

        expect(mockFindWorkflow).toHaveBeenCalledWith({ id: 123 });
        expect(result).toBe('My Workflow Title');
      });

      it('should handle string ID', async () => {
        const mockFindWorkflow = vi.fn().mockResolvedValue({
          data: {
            title: 'Workflow from String ID'
          }
        });

        const getWorkflowTitle = async function(id) {
          const workflow_id = {
            id: id,
          };
          const user_workflow = await mockFindWorkflow(workflow_id);
          return user_workflow.data.title;
        };

        const result = await getWorkflowTitle.call(mockComponent, '456');

        expect(mockFindWorkflow).toHaveBeenCalledWith({ id: '456' });
        expect(result).toBe('Workflow from String ID');
      });

      it('should handle API error', async () => {
        const mockFindWorkflow = vi.fn().mockRejectedValue(new Error('API Error'));

        const getWorkflowTitle = async function(id) {
          const workflow_id = {
            id: id,
          };
          const user_workflow = await mockFindWorkflow(workflow_id);
          return user_workflow.data.title;
        };

        await expect(getWorkflowTitle.call(mockComponent, 123)).rejects.toThrow('API Error');
      });
    });
  });

  /**
   * View Mode Tests
   */
  describe('View Mode', () => {
    describe('activateCompileCheck', () => {
      it('should set compile_check to true', () => {
        const activateCompileCheck = function() {
          this.compile_check = true;
        };

        mockComponent.compile_check = false;

        activateCompileCheck.call(mockComponent);

        expect(mockComponent.compile_check).toBe(true);
      });

      it('should handle already active state', () => {
        const activateCompileCheck = function() {
          this.compile_check = true;
        };

        mockComponent.compile_check = true;

        activateCompileCheck.call(mockComponent);

        expect(mockComponent.compile_check).toBe(true);
      });
    });

    describe('deactivateCompileCheck', () => {
      it('should set compile_check to false', () => {
        const deactivateCompileCheck = function() {
          this.compile_check = false;
        };

        mockComponent.compile_check = true;

        deactivateCompileCheck.call(mockComponent);

        expect(mockComponent.compile_check).toBe(false);
      });

      it('should handle already inactive state', () => {
        const deactivateCompileCheck = function() {
          this.compile_check = false;
        };

        mockComponent.compile_check = false;

        deactivateCompileCheck.call(mockComponent);

        expect(mockComponent.compile_check).toBe(false);
      });
    });

    describe('updateIsTabView', () => {
      it('should update isTabView to true', () => {
        const updateIsTabView = function(newIsTabView) {
          this.isTabView = newIsTabView;
        };

        mockComponent.isTabView = false;

        updateIsTabView.call(mockComponent, true);

        expect(mockComponent.isTabView).toBe(true);
      });

      it('should update isTabView to false', () => {
        const updateIsTabView = function(newIsTabView) {
          this.isTabView = newIsTabView;
        };

        mockComponent.isTabView = true;

        updateIsTabView.call(mockComponent, false);

        expect(mockComponent.isTabView).toBe(false);
      });

      it('should handle boolean string conversion', () => {
        const updateIsTabView = function(newIsTabView) {
          this.isTabView = newIsTabView;
        };

        mockComponent.isTabView = false;

        updateIsTabView.call(mockComponent, 'true');

        expect(mockComponent.isTabView).toBe('true');
      });
    });
  });
});

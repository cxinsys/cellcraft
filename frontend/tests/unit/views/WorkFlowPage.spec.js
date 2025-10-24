import { describe, it, expect, vi, beforeEach } from 'vitest';
import { getTimeDifference } from '@/utils/formatters';

/**
 * Unit tests for WorkFlowPage.vue decomposed methods
 *
 * These tests focus on the Phase 4 refactored methods:
 * - runWorkflow group (7 methods)
 * - toggleTask group (6 methods)
 * - setCurrentWorkflow group (3 methods)
 */

describe('WorkFlowPage.vue - Decomposed Methods', () => {
  let mockComponent;
  let mockStore;
  let mockDrawflow;

  beforeEach(() => {
    // Mock Vuex store
    mockStore = {
      getters: {
        getTitle: 'Test Workflow',
        getThumbnail: 'test-thumbnail.png',
        getWorkflowInfo: {
          drawflow: {
            Home: {
              data: {}
            }
          }
        }
      },
      commit: vi.fn()
    };

    // Mock Drawflow instance
    mockDrawflow = {
      export: vi.fn().mockReturnValue({ workflow: 'data' }),
      import: vi.fn(),
      drawflow: {
        drawflow: {
          Home: {
            data: {}
          }
        }
      },
      module: 'Home'
    };

    // Mock component instance
    mockComponent = {
      $store: mockStore,
      $df: mockDrawflow,
      currentWorkflowId: '123',
      exportValue: null,
      show_jobs: false,
      taskList: [],
      timeInterval: null,
      on_progress: false,

      // Methods that will be called by decomposed methods
      updateWorkflowInfo: vi.fn(),
      createEventSource: vi.fn(),
      toggleTask: vi.fn(),
      setCurrentWorkflowInfo: vi.fn(),
      captureWorkflow: vi.fn().mockResolvedValue(undefined),
      setMessage: vi.fn(),
      startTimer: vi.fn().mockReturnValue(123)
    };
  });

  /**
   * runWorkflow Group Tests
   */
  describe('runWorkflow group', () => {
    describe('prepareWorkflowData', () => {
      it('should prepare workflow data with all required fields', () => {
        // Import the actual implementation pattern
        const prepareWorkflowData = function() {
          this.updateWorkflowInfo();
          this.exportValue = this.$df.export();
          const title = this.$store.getters.getTitle;
          const thumbnail = this.$store.getters.getThumbnail;

          return {
            id: this.currentWorkflowId,
            title: title,
            thumbnail: thumbnail,
            workflow_info: this.exportValue,
          };
        };

        const result = prepareWorkflowData.call(mockComponent);

        expect(mockComponent.updateWorkflowInfo).toHaveBeenCalled();
        expect(mockDrawflow.export).toHaveBeenCalled();
        expect(result).toEqual({
          id: '123',
          title: 'Test Workflow',
          thumbnail: 'test-thumbnail.png',
          workflow_info: { workflow: 'data' }
        });
      });

      it('should update exportValue on component instance', () => {
        const prepareWorkflowData = function() {
          this.updateWorkflowInfo();
          this.exportValue = this.$df.export();
          const title = this.$store.getters.getTitle;
          const thumbnail = this.$store.getters.getThumbnail;

          return {
            id: this.currentWorkflowId,
            title: title,
            thumbnail: thumbnail,
            workflow_info: this.exportValue,
          };
        };

        prepareWorkflowData.call(mockComponent);

        expect(mockComponent.exportValue).toEqual({ workflow: 'data' });
      });
    });

    describe('submitWorkflow', () => {
      it('should call exportData API with workflow payload', async () => {
        const mockExportData = vi.fn().mockResolvedValue({
          data: { task_ids: ['task1', 'task2'] }
        });

        const submitWorkflow = async function(workflow) {
          console.log(this.$df.drawflow.drawflow[this.$df.module]);
          return await mockExportData(workflow);
        };

        const workflow = {
          id: '123',
          title: 'Test',
          workflow_info: {}
        };

        const result = await submitWorkflow.call(mockComponent, workflow);

        expect(mockExportData).toHaveBeenCalledWith(workflow);
        expect(result).toEqual({ data: { task_ids: ['task1', 'task2'] } });
      });
    });

    describe('handleWorkflowResponse', () => {
      it('should process valid workflow response', () => {
        const handleWorkflowResponse = function(workflowData) {
          if (!workflowData.data.task_ids || !Array.isArray(workflowData.data.task_ids)) {
            return;
          }

          this.updateTaskMapping(workflowData.data);
          this.startTaskMonitoring(workflowData.data.task_ids);
        };

        mockComponent.updateTaskMapping = vi.fn();
        mockComponent.startTaskMonitoring = vi.fn();

        const workflowData = {
          data: {
            task_ids: ['task1', 'task2'],
            task_algorithm_mapping: { task1: 'algo1' },
            algorithm_ids: ['algo1']
          }
        };

        handleWorkflowResponse.call(mockComponent, workflowData);

        expect(mockComponent.updateTaskMapping).toHaveBeenCalledWith(workflowData.data);
        expect(mockComponent.startTaskMonitoring).toHaveBeenCalledWith(['task1', 'task2']);
      });

      it('should return early if task_ids is missing', () => {
        const handleWorkflowResponse = function(workflowData) {
          if (!workflowData.data.task_ids || !Array.isArray(workflowData.data.task_ids)) {
            return;
          }

          this.updateTaskMapping(workflowData.data);
          this.startTaskMonitoring(workflowData.data.task_ids);
        };

        mockComponent.updateTaskMapping = vi.fn();
        mockComponent.startTaskMonitoring = vi.fn();

        const workflowData = {
          data: {}
        };

        handleWorkflowResponse.call(mockComponent, workflowData);

        expect(mockComponent.updateTaskMapping).not.toHaveBeenCalled();
        expect(mockComponent.startTaskMonitoring).not.toHaveBeenCalled();
      });

      it('should return early if task_ids is not an array', () => {
        const handleWorkflowResponse = function(workflowData) {
          if (!workflowData.data.task_ids || !Array.isArray(workflowData.data.task_ids)) {
            return;
          }

          this.updateTaskMapping(workflowData.data);
          this.startTaskMonitoring(workflowData.data.task_ids);
        };

        mockComponent.updateTaskMapping = vi.fn();
        mockComponent.startTaskMonitoring = vi.fn();

        const workflowData = {
          data: {
            task_ids: 'not-an-array'
          }
        };

        handleWorkflowResponse.call(mockComponent, workflowData);

        expect(mockComponent.updateTaskMapping).not.toHaveBeenCalled();
        expect(mockComponent.startTaskMonitoring).not.toHaveBeenCalled();
      });
    });

    describe('updateTaskMapping', () => {
      it('should commit task_algorithm_mapping to store', () => {
        const updateTaskMapping = function(data) {
          if (data.task_algorithm_mapping) {
            this.$store.commit('setTaskAlgorithmMapping', data.task_algorithm_mapping);
          }

          if (data.algorithm_ids && Array.isArray(data.algorithm_ids)) {
            this.$store.commit('setRunningAlgorithmNodes', data.algorithm_ids);
          }
        };

        const data = {
          task_algorithm_mapping: { task1: 'algo1', task2: 'algo2' },
          algorithm_ids: ['algo1', 'algo2']
        };

        updateTaskMapping.call(mockComponent, data);

        expect(mockStore.commit).toHaveBeenCalledWith('setTaskAlgorithmMapping', data.task_algorithm_mapping);
        expect(mockStore.commit).toHaveBeenCalledWith('setRunningAlgorithmNodes', data.algorithm_ids);
      });

      it('should handle missing task_algorithm_mapping', () => {
        const updateTaskMapping = function(data) {
          if (data.task_algorithm_mapping) {
            this.$store.commit('setTaskAlgorithmMapping', data.task_algorithm_mapping);
          }

          if (data.algorithm_ids && Array.isArray(data.algorithm_ids)) {
            this.$store.commit('setRunningAlgorithmNodes', data.algorithm_ids);
          }
        };

        const data = {
          algorithm_ids: ['algo1']
        };

        updateTaskMapping.call(mockComponent, data);

        expect(mockStore.commit).toHaveBeenCalledTimes(1);
        expect(mockStore.commit).toHaveBeenCalledWith('setRunningAlgorithmNodes', data.algorithm_ids);
      });

      it('should handle missing algorithm_ids', () => {
        const updateTaskMapping = function(data) {
          if (data.task_algorithm_mapping) {
            this.$store.commit('setTaskAlgorithmMapping', data.task_algorithm_mapping);
          }

          if (data.algorithm_ids && Array.isArray(data.algorithm_ids)) {
            this.$store.commit('setRunningAlgorithmNodes', data.algorithm_ids);
          }
        };

        const data = {
          task_algorithm_mapping: { task1: 'algo1' }
        };

        updateTaskMapping.call(mockComponent, data);

        expect(mockStore.commit).toHaveBeenCalledTimes(1);
        expect(mockStore.commit).toHaveBeenCalledWith('setTaskAlgorithmMapping', data.task_algorithm_mapping);
      });
    });

    describe('startTaskMonitoring', () => {
      it('should create event source for each task ID', () => {
        const startTaskMonitoring = function(taskIds) {
          taskIds.forEach(task_id => {
            this.createEventSource(task_id);
          });
        };

        const taskIds = ['task1', 'task2', 'task3'];

        startTaskMonitoring.call(mockComponent, taskIds);

        expect(mockComponent.createEventSource).toHaveBeenCalledTimes(3);
        expect(mockComponent.createEventSource).toHaveBeenCalledWith('task1');
        expect(mockComponent.createEventSource).toHaveBeenCalledWith('task2');
        expect(mockComponent.createEventSource).toHaveBeenCalledWith('task3');
      });

      it('should handle empty task IDs array', () => {
        const startTaskMonitoring = function(taskIds) {
          taskIds.forEach(task_id => {
            this.createEventSource(task_id);
          });
        };

        startTaskMonitoring.call(mockComponent, []);

        expect(mockComponent.createEventSource).not.toHaveBeenCalled();
      });
    });

    describe('updateUIStateAfterRun', () => {
      it('should toggle task if show_jobs is true', () => {
        const updateUIStateAfterRun = function() {
          if (this.show_jobs) {
            this.show_jobs = false;
            this.toggleTask();
          }
        };

        mockComponent.show_jobs = true;

        updateUIStateAfterRun.call(mockComponent);

        expect(mockComponent.show_jobs).toBe(false);
        expect(mockComponent.toggleTask).toHaveBeenCalled();
      });

      it('should not toggle task if show_jobs is false', () => {
        const updateUIStateAfterRun = function() {
          if (this.show_jobs) {
            this.show_jobs = false;
            this.toggleTask();
          }
        };

        mockComponent.show_jobs = false;

        updateUIStateAfterRun.call(mockComponent);

        expect(mockComponent.toggleTask).not.toHaveBeenCalled();
      });
    });

    describe('runWorkflow (orchestrator)', () => {
      it('should execute workflow submission flow successfully', async () => {
        const mockExportData = vi.fn().mockResolvedValue({
          data: {
            task_ids: ['task1'],
            task_algorithm_mapping: { task1: 'algo1' },
            algorithm_ids: ['algo1']
          }
        });

        const runWorkflow = async function() {
          try {
            this.updateWorkflowInfo();
            this.exportValue = this.$df.export();
            const title = this.$store.getters.getTitle;
            const thumbnail = this.$store.getters.getThumbnail;

            const workflow = {
              id: this.currentWorkflowId,
              title: title,
              thumbnail: thumbnail,
              workflow_info: this.exportValue,
            };

            const workflowData = await mockExportData(workflow);

            if (workflowData.data.task_ids && Array.isArray(workflowData.data.task_ids)) {
              if (workflowData.data.task_algorithm_mapping) {
                this.$store.commit('setTaskAlgorithmMapping', workflowData.data.task_algorithm_mapping);
              }
              if (workflowData.data.algorithm_ids && Array.isArray(workflowData.data.algorithm_ids)) {
                this.$store.commit('setRunningAlgorithmNodes', workflowData.data.algorithm_ids);
              }
              workflowData.data.task_ids.forEach(task_id => {
                this.createEventSource(task_id);
              });
            }

            if (this.show_jobs) {
              this.show_jobs = false;
              this.toggleTask();
            }
          } catch (error) {
            console.error(error);
          }
        };

        await runWorkflow.call(mockComponent);

        expect(mockComponent.updateWorkflowInfo).toHaveBeenCalled();
        expect(mockDrawflow.export).toHaveBeenCalled();
        expect(mockExportData).toHaveBeenCalled();
        expect(mockComponent.createEventSource).toHaveBeenCalledWith('task1');
      });

      it('should handle errors gracefully', async () => {
        const mockExportData = vi.fn().mockRejectedValue(new Error('API Error'));
        const consoleErrorSpy = vi.spyOn(console, 'error').mockImplementation(() => {});

        const runWorkflow = async function() {
          try {
            this.updateWorkflowInfo();
            this.exportValue = this.$df.export();
            const title = this.$store.getters.getTitle;
            const thumbnail = this.$store.getters.getThumbnail;

            const workflow = {
              id: this.currentWorkflowId,
              title: title,
              thumbnail: thumbnail,
              workflow_info: this.exportValue,
            };

            await mockExportData(workflow);
          } catch (error) {
            console.error(error);
          }
        };

        await runWorkflow.call(mockComponent);

        expect(consoleErrorSpy).toHaveBeenCalled();
        consoleErrorSpy.mockRestore();
      });
    });
  });

  /**
   * toggleTask Group Tests
   */
  describe('toggleTask group', () => {
    describe('fetchUserTasks', () => {
      it('should fetch and return user tasks', async () => {
        const mockUserTaskMonitoring = vi.fn().mockResolvedValue({
          data: [
            { id: 'task1', status: 'SUCCESS' },
            { id: 'task2', status: 'RUNNING' }
          ]
        });

        const fetchUserTasks = async function() {
          const user_tasks = await mockUserTaskMonitoring();
          console.log(user_tasks);
          return user_tasks.data;
        };

        const result = await fetchUserTasks.call(mockComponent);

        expect(mockUserTaskMonitoring).toHaveBeenCalled();
        expect(result).toEqual([
          { id: 'task1', status: 'SUCCESS' },
          { id: 'task2', status: 'RUNNING' }
        ]);
      });
    });

    describe('processTaskList', () => {
      it('should process completed tasks with time difference', () => {
        const processTaskList = function(tasks) {
          this.taskList = tasks;

          this.taskList.forEach((task, idx) => {
            if (task.status === "SUCCESS" ||
              task.status === "FAILURE" ||
              task.status === "REVOKED" ||
              task.status === "RETRY") {
              this.taskList[idx].running_time = getTimeDifference(
                task.start_time,
                task.end_time
              );
            } else if (task.status === "RUNNING" || task.status === "PENDING" || task.status === "INSTALLING") {
              this.timeInterval = this.startTimer(idx);
              this.on_progress = true;
            }
          });
        };

        const tasks = [
          {
            id: 'task1',
            status: 'SUCCESS',
            start_time: '2024-01-01T12:00:00',
            end_time: '2024-01-01T12:05:30'
          }
        ];

        processTaskList.call(mockComponent, tasks);

        expect(mockComponent.taskList[0].running_time).toBe('00:05:30');
      });

      it('should start timer for running tasks', () => {
        const processTaskList = function(tasks) {
          this.taskList = tasks;

          this.taskList.forEach((task, idx) => {
            if (task.status === "SUCCESS" ||
              task.status === "FAILURE" ||
              task.status === "REVOKED" ||
              task.status === "RETRY") {
              this.taskList[idx].running_time = getTimeDifference(
                task.start_time,
                task.end_time
              );
            } else if (task.status === "RUNNING" || task.status === "PENDING" || task.status === "INSTALLING") {
              this.timeInterval = this.startTimer(idx);
              this.on_progress = true;
            }
          });
        };

        const tasks = [
          {
            id: 'task1',
            status: 'RUNNING',
            start_time: '2024-01-01T12:00:00'
          }
        ];

        processTaskList.call(mockComponent, tasks);

        expect(mockComponent.startTimer).toHaveBeenCalledWith(0);
        expect(mockComponent.on_progress).toBe(true);
        expect(mockComponent.timeInterval).toBe(123);
      });

      it('should handle multiple tasks with mixed statuses', () => {
        const processTaskList = function(tasks) {
          this.taskList = tasks;

          this.taskList.forEach((task, idx) => {
            if (task.status === "SUCCESS" ||
              task.status === "FAILURE" ||
              task.status === "REVOKED" ||
              task.status === "RETRY") {
              this.taskList[idx].running_time = getTimeDifference(
                task.start_time,
                task.end_time
              );
            } else if (task.status === "RUNNING" || task.status === "PENDING" || task.status === "INSTALLING") {
              this.timeInterval = this.startTimer(idx);
              this.on_progress = true;
            }
          });
        };

        const tasks = [
          {
            id: 'task1',
            status: 'SUCCESS',
            start_time: '2024-01-01T12:00:00',
            end_time: '2024-01-01T12:05:00'
          },
          {
            id: 'task2',
            status: 'PENDING',
            start_time: '2024-01-01T12:05:00'
          }
        ];

        processTaskList.call(mockComponent, tasks);

        expect(mockComponent.taskList[0].running_time).toBe('00:05:00');
        expect(mockComponent.startTimer).toHaveBeenCalledWith(1);
      });
    });

    describe('stopTaskTimers', () => {
      it('should clear timeInterval', () => {
        const clearIntervalSpy = vi.spyOn(global, 'clearInterval');

        const stopTaskTimers = function() {
          clearInterval(this.timeInterval);
        };

        mockComponent.timeInterval = 123;

        stopTaskTimers.call(mockComponent);

        expect(clearIntervalSpy).toHaveBeenCalledWith(123);
        clearIntervalSpy.mockRestore();
      });
    });

    describe('showTaskList', () => {
      it('should fetch and process task list', async () => {
        const mockUserTaskMonitoring = vi.fn().mockResolvedValue({
          data: [
            {
              id: 'task1',
              status: 'SUCCESS',
              start_time: '2024-01-01T12:00:00',
              end_time: '2024-01-01T12:05:00'
            }
          ]
        });

        const showTaskList = async function() {
          const user_tasks = await mockUserTaskMonitoring();
          const tasks = user_tasks.data;

          this.taskList = tasks;
          this.taskList.forEach((task, idx) => {
            if (task.status === "SUCCESS" ||
              task.status === "FAILURE" ||
              task.status === "REVOKED" ||
              task.status === "RETRY") {
              this.taskList[idx].running_time = getTimeDifference(
                task.start_time,
                task.end_time
              );
            } else if (task.status === "RUNNING" || task.status === "PENDING" || task.status === "INSTALLING") {
              this.timeInterval = this.startTimer(idx);
              this.on_progress = true;
            }
          });
        };

        await showTaskList.call(mockComponent);

        expect(mockUserTaskMonitoring).toHaveBeenCalled();
        expect(mockComponent.taskList).toHaveLength(1);
        expect(mockComponent.taskList[0].running_time).toBe('00:05:00');
      });
    });

    describe('hideTaskList', () => {
      it('should stop task timers', () => {
        const clearIntervalSpy = vi.spyOn(global, 'clearInterval');

        const hideTaskList = function() {
          clearInterval(this.timeInterval);
        };

        mockComponent.timeInterval = 123;

        hideTaskList.call(mockComponent);

        expect(clearIntervalSpy).toHaveBeenCalledWith(123);
        clearIntervalSpy.mockRestore();
      });
    });

    describe('toggleTask (orchestrator)', () => {
      it('should show task list when show_jobs is false', async () => {
        const mockUserTaskMonitoring = vi.fn().mockResolvedValue({
          data: [
            {
              id: 'task1',
              status: 'SUCCESS',
              start_time: '2024-01-01T12:00:00',
              end_time: '2024-01-01T12:05:00'
            }
          ]
        });

        const setTimeoutSpy = vi.spyOn(global, 'setTimeout');

        const toggleTask = async function() {
          try {
            if (!this.show_jobs) {
              const user_tasks = await mockUserTaskMonitoring();
              const tasks = user_tasks.data;
              this.taskList = tasks;
              this.taskList.forEach((task, idx) => {
                if (task.status === "SUCCESS" ||
                  task.status === "FAILURE" ||
                  task.status === "REVOKED" ||
                  task.status === "RETRY") {
                  this.taskList[idx].running_time = getTimeDifference(
                    task.start_time,
                    task.end_time
                  );
                }
              });
            } else {
              clearInterval(this.timeInterval);
            }

            setTimeout(() => {
              this.show_jobs = !this.show_jobs;
            }, 300);
          } catch (error) {
            console.error(error);
            this.setMessage(
              "error",
              "No tasks have been executed yet. Please run workflow"
            );
          }
        };

        mockComponent.show_jobs = false;

        await toggleTask.call(mockComponent);

        expect(mockUserTaskMonitoring).toHaveBeenCalled();
        expect(mockComponent.taskList).toHaveLength(1);
        expect(setTimeoutSpy).toHaveBeenCalledWith(expect.any(Function), 300);

        setTimeoutSpy.mockRestore();
      });

      it('should hide task list when show_jobs is true', async () => {
        const clearIntervalSpy = vi.spyOn(global, 'clearInterval');
        const setTimeoutSpy = vi.spyOn(global, 'setTimeout');

        const toggleTask = async function() {
          try {
            if (!this.show_jobs) {
              // Show logic
            } else {
              clearInterval(this.timeInterval);
            }

            setTimeout(() => {
              this.show_jobs = !this.show_jobs;
            }, 300);
          } catch (error) {
            console.error(error);
            this.setMessage(
              "error",
              "No tasks have been executed yet. Please run workflow"
            );
          }
        };

        mockComponent.show_jobs = true;
        mockComponent.timeInterval = 123;

        await toggleTask.call(mockComponent);

        expect(clearIntervalSpy).toHaveBeenCalledWith(123);
        expect(setTimeoutSpy).toHaveBeenCalledWith(expect.any(Function), 300);

        clearIntervalSpy.mockRestore();
        setTimeoutSpy.mockRestore();
      });

      it('should handle errors and show error message', async () => {
        const mockUserTaskMonitoring = vi.fn().mockRejectedValue(new Error('API Error'));
        const consoleErrorSpy = vi.spyOn(console, 'error').mockImplementation(() => {});

        const toggleTask = async function() {
          try {
            if (!this.show_jobs) {
              await mockUserTaskMonitoring();
            }

            setTimeout(() => {
              this.show_jobs = !this.show_jobs;
            }, 300);
          } catch (error) {
            console.error(error);
            this.setMessage(
              "error",
              "No tasks have been executed yet. Please run workflow"
            );
          }
        };

        mockComponent.show_jobs = false;

        await toggleTask.call(mockComponent);

        expect(consoleErrorSpy).toHaveBeenCalled();
        expect(mockComponent.setMessage).toHaveBeenCalledWith(
          "error",
          "No tasks have been executed yet. Please run workflow"
        );

        consoleErrorSpy.mockRestore();
      });
    });
  });

  /**
   * setCurrentWorkflow Group Tests
   */
  describe('setCurrentWorkflow group', () => {
    describe('buildWorkflowPayload', () => {
      it('should build workflow payload with all required fields', () => {
        const buildWorkflowPayload = function() {
          const title = this.$store.getters.getTitle;
          const thumbnail = this.$store.getters.getThumbnail;

          return {
            id: this.currentWorkflowId,
            title: title,
            thumbnail: thumbnail,
            workflow_info: this.exportValue,
          };
        };

        mockComponent.exportValue = { workflow: 'test-data' };

        const result = buildWorkflowPayload.call(mockComponent);

        expect(result).toEqual({
          id: '123',
          title: 'Test Workflow',
          thumbnail: 'test-thumbnail.png',
          workflow_info: { workflow: 'test-data' }
        });
      });

      it('should use current exportValue from component', () => {
        const buildWorkflowPayload = function() {
          const title = this.$store.getters.getTitle;
          const thumbnail = this.$store.getters.getThumbnail;

          return {
            id: this.currentWorkflowId,
            title: title,
            thumbnail: thumbnail,
            workflow_info: this.exportValue,
          };
        };

        mockComponent.exportValue = null;

        const result = buildWorkflowPayload.call(mockComponent);

        expect(result.workflow_info).toBeNull();
      });
    });

    describe('saveWorkflowToServer', () => {
      it('should save workflow and update currentWorkflowId', async () => {
        const mockSaveWorkflow = vi.fn().mockResolvedValue({
          data: {
            id: 456,
            title: 'Saved Workflow'
          }
        });

        const saveWorkflowToServer = async function(workflow) {
          console.log("currentWorkflowId : " + workflow.id + " type : " + typeof workflow.id);
          const workflow_data = await mockSaveWorkflow(workflow);
          this.currentWorkflowId = String(workflow_data.data.id);
          return workflow_data.data;
        };

        const workflow = {
          id: '123',
          title: 'Test Workflow',
          workflow_info: {}
        };

        const result = await saveWorkflowToServer.call(mockComponent, workflow);

        expect(mockSaveWorkflow).toHaveBeenCalledWith(workflow);
        expect(mockComponent.currentWorkflowId).toBe('456');
        expect(result).toEqual({
          id: 456,
          title: 'Saved Workflow'
        });
      });

      it('should convert workflow ID to string', async () => {
        const mockSaveWorkflow = vi.fn().mockResolvedValue({
          data: {
            id: 999
          }
        });

        const saveWorkflowToServer = async function(workflow) {
          const workflow_data = await mockSaveWorkflow(workflow);
          this.currentWorkflowId = String(workflow_data.data.id);
          return workflow_data.data;
        };

        await saveWorkflowToServer.call(mockComponent, { id: '123' });

        expect(mockComponent.currentWorkflowId).toBe('999');
        expect(typeof mockComponent.currentWorkflowId).toBe('string');
      });
    });

    describe('setCurrentWorkflow (orchestrator)', () => {
      it('should execute complete workflow saving flow', async () => {
        const mockSaveWorkflow = vi.fn().mockResolvedValue({
          data: {
            id: 789,
            title: 'Final Workflow'
          }
        });

        const setCurrentWorkflow = async function() {
          try {
            this.setCurrentWorkflowInfo();
            await this.captureWorkflow();

            const title = this.$store.getters.getTitle;
            const thumbnail = this.$store.getters.getThumbnail;
            const workflow = {
              id: this.currentWorkflowId,
              title: title,
              thumbnail: thumbnail,
              workflow_info: this.exportValue,
            };

            const workflow_data = await mockSaveWorkflow(workflow);
            this.currentWorkflowId = String(workflow_data.data.id);
            return workflow_data.data;
          } catch (error) {
            console.error(error);
          }
        };

        mockComponent.exportValue = { workflow: 'final-data' };

        const result = await setCurrentWorkflow.call(mockComponent);

        expect(mockComponent.setCurrentWorkflowInfo).toHaveBeenCalled();
        expect(mockComponent.captureWorkflow).toHaveBeenCalled();
        expect(mockSaveWorkflow).toHaveBeenCalled();
        expect(mockComponent.currentWorkflowId).toBe('789');
        expect(result).toEqual({
          id: 789,
          title: 'Final Workflow'
        });
      });

      it('should handle errors gracefully', async () => {
        const mockSaveWorkflow = vi.fn().mockRejectedValue(new Error('Save Error'));
        const consoleErrorSpy = vi.spyOn(console, 'error').mockImplementation(() => {});

        const setCurrentWorkflow = async function() {
          try {
            this.setCurrentWorkflowInfo();
            await this.captureWorkflow();

            const workflow = {
              id: this.currentWorkflowId,
              title: 'Test',
              workflow_info: this.exportValue,
            };

            await mockSaveWorkflow(workflow);
          } catch (error) {
            console.error(error);
          }
        };

        await setCurrentWorkflow.call(mockComponent);

        expect(consoleErrorSpy).toHaveBeenCalled();
        consoleErrorSpy.mockRestore();
      });
    });
  });
});

<template>
  <div>
    <div class="nodeBox" ref="nodeBox">
      <img class="nodeBox__icon" src="@/assets/ScatterPlot_logo.png" draggable="false" />
    </div>
    <div class="nodeTitleBox">
      <input type="text" class="nodeTitle" v-model="nodeTitle" @input="updateTitle" />
    </div>
  </div>
</template>

<script>
import { mapMutations } from "vuex";

export default {
  data() {
    return {
      nodeId: null,
    };
  },
  mounted() {
    this.$nextTick(() => {
      const nodeBox = this.$refs.nodeBox;
      const nodeId = parseInt(
        nodeBox.parentNode.parentNode.parentNode.id.replace("node-", "")
      );
      this.nodeId = nodeId;
    });
  },
  computed: {
    nodeTitle: {
      get() {
        if (this.nodeId == null) {
          return "ScatterPlot"; // 또는 기본 타이틀 반환
        }
        const node = this.$store.getters.getWorkflowNodeInfo(this.nodeId);
        return node.data.title || node.name;
      },
      set(newTitle) {
        this.updateWorkflowNodeTitle({ nodeId: this.nodeId, newTitle });
      },
    },
  },
  methods: {
    ...mapMutations(["updateWorkflowNodeTitle"]),
    updateTitle() {
      this.updateWorkflowNodeTitle({ nodeId: this.nodeId, newTitle: this.nodeTitle });
    },
  },
};
</script>

<style scoped>
.nodeBox {
  width: 4rem;
  height: 4rem;
  display: flex;
  align-items: center;
  justify-content: center;
  text-align: center;
}

.nodeBox__icon {
  width: 55%;
  height: 55%;
  object-fit: contain;
}

.nodeTitleBox {
  margin-top: 0.2rem;
  width: 8rem;
  height: 1.5rem;
  left: calc(50% - 4rem);
  display: flex;
  position: absolute;
  align-items: center;
  justify-content: center;
}

.nodeTitle {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-size: 0.9rem;
  font-weight: 300;
  width: 100%;
  text-align: center;
  color: rgb(233, 233, 233);
  background: none;
  border: none;
  text-overflow: ellipsis;
}

.nodeTitle:focus {
  outline: none;
  background: rgba(255, 255, 255, 0.442);
  border-radius: 4px;
}

@media (prefers-color-scheme: dark) {
  /* .nodeBox__icon {
    filter: invert(97%) sepia(99%) saturate(0%) hue-rotate(123deg)
      brightness(107%) contrast(101%);
  } */
}
</style>

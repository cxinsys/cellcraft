<template>
  <div>
    <div class="nodeBox" ref="nodeBox">
      <img class="nodeBox__icon" src="@/assets/input_file.png" draggable="false" />
    </div>
    <div class="nodeTitleBox">
      <input type="text" class="nodeTitle" v-model="nodeTitle" @input="updateTitle" df-title />
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
          return "Input File"; // 또는 기본 타이틀 반환
        }
        const node = this.$store.getters.getNodeInfo(this.nodeId);
        return node.title;
      },
      set(newTitle) {
        this.updateNodeTitle({ nodeId: this.nodeId, newTitle });
      },
    },
  },
  methods: {
    ...mapMutations(["updateNodeTitle"]),
    updateTitle() {
      this.updateNodeTitle({ nodeId: this.nodeId, newTitle: this.nodeTitle });
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

.layout {
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
}

.toggle-layout {
  margin: 20px;
  border-radius: 10px;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
}

.toggle-layout__row {
  width: 100%;
  height: auto;
  display: flex;
  align-items: center;
  justify-content: space-between;
}

.fileUpload-form {
  width: 100%;
  height: 100%;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  text-align: center;
}

.fileUpload-form__button {
  display: inline-block;
  width: 46%;
  height: 1.5rem;
  padding-top: 10px;
  border: 1px solid black;
  color: black;
  background: white;
  margin-bottom: 5px;
  cursor: pointer;
  vertical-align: middle;
  font-size: 14px;
}

.fileUpload-form input[type="file"],
.fileUpload-form input[type="submit"] {
  position: absolute;
  width: 0;
  height: 0;
  padding: 0;
  overflow: hidden;
  border: 0;
}

.fileUpload-form__title {
  font-size: 1rem;
  font-weight: bold;
}

/* .fileUpload-form__info{
  
  } */

.fileUpload-form__arrow {
  color: black;
  width: 1.5rem;
  height: auto;
}

.toggle-icon {
  background: rgb(227, 243, 252);
  width: 20px;
  height: 20px;
  border-radius: 50%;
  margin-top: 5px;
  display: flex;
  align-items: center;
  justify-content: center;
  text-align: center;
}

.nodeTitleBox {
  margin-top: 0.2rem;
  width: 10rem;
  height: 1.5rem;
  left: calc(50% - 5rem);
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